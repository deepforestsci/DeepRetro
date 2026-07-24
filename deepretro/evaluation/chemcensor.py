"""Optional adapter for the official local ChemCensor implementation.

DeepRetro intentionally delegates precedent matching, atom mapping, reaction
centre/context extraction, and functional-group compatibility to the official
engine.  The optional dependency is imported only when an operator explicitly
constructs :class:`OfficialChemCensorScorer`.
"""

from __future__ import annotations

import importlib
import logging
import sys
from importlib import metadata
from pathlib import Path
from types import ModuleType
from typing import Any, Callable

from deepretro.evaluation.cache import ScoreCache
from deepretro.evaluation.io import file_sha256
from deepretro.evaluation.types import (
    CanonicalPrediction,
    ChemCensorProvenance,
    ChemCensorScore,
)

logger = logging.getLogger(__name__)


class ChemCensorUnavailableError(RuntimeError):
    """Raised when the optional official scorer cannot be initialized."""


class OfficialChemCensorScorer:
    """Score canonical DeepRetro predictions with a local official engine.

    Parameters are deliberately a subset of the official constructor.  Tests
    inject a module loader; production callers supply an approved local SQLite
    database and a compatible ``chemcensor`` installation.
    """

    def __init__(
        self,
        database_path: str | Path,
        *,
        database_version: str,
        max_center_type: int = 4,
        find_exact_match: bool = True,
        cache_path: str | Path | None = None,
        module_loader: Callable[[str], ModuleType] = importlib.import_module,
        package_version_getter: Callable[[str], str] = metadata.version,
    ) -> None:
        """Initialize the official scorer and immutable provenance payload."""
        if sys.version_info < (3, 12):
            raise ChemCensorUnavailableError(
                "the optional ChemCensor package requires Python 3.12 or later"
            )
        if isinstance(max_center_type, bool) or not isinstance(max_center_type, int):
            raise ValueError("max_center_type must be an integer")
        if max_center_type < 1:
            raise ValueError("max_center_type must be positive")
        if not isinstance(find_exact_match, bool):
            raise ValueError("find_exact_match must be a boolean")

        self._database_path = Path(database_path).expanduser().resolve()
        if not self._database_path.is_file():
            raise FileNotFoundError(
                f"ChemCensor database does not exist: {self._database_path}"
            )
        try:
            module = module_loader("chemcensor")
        except ModuleNotFoundError as exc:
            dependency_hint = (
                "the optional ChemCensor package is not installed"
                if exc.name in {None, "chemcensor"}
                else f"the optional ChemCensor dependency {exc.name!r} is missing"
            )
            raise ChemCensorUnavailableError(
                f"{dependency_hint}; install a compatible local official release "
                "to use score"
            ) from exc

        factory = getattr(module, "ChemCensor", None)
        if not callable(factory):
            raise ChemCensorUnavailableError(
                "the installed ChemCensor package does not expose ChemCensor"
            )
        try:
            self._engine = factory(
                db_path=self._database_path,
                max_center_type=max_center_type,
                find_exact_match=find_exact_match,
            )
        except Exception as exc:  # external database/configuration boundary
            raise ChemCensorUnavailableError(
                f"failed to initialize the official ChemCensor engine: {exc}"
            ) from exc
        if not callable(getattr(self._engine, "evaluate", None)):
            raise ChemCensorUnavailableError(
                "the installed ChemCensor release does not provide evaluate(); "
                "a dual-score release is required"
            )

        try:
            package_version = package_version_getter("chemcensor")
        except metadata.PackageNotFoundError:
            package_version = "unknown"
        self.provenance = ChemCensorProvenance(
            scorer_name="chemcensor",
            package_version=package_version,
            database_filename=self._database_path.name,
            database_version=database_version,
            database_sha256=file_sha256(self._database_path),
            database_size_bytes=self._database_path.stat().st_size,
            max_center_type=max_center_type,
            find_exact_match=find_exact_match,
        )
        self._cache: ScoreCache | None = None
        if cache_path is not None:
            try:
                self._cache = ScoreCache(cache_path)
            except Exception as exc:  # cache is a performance optimization only
                logger.warning(
                    "ChemCensor score cache initialization failed; bypassing cache: %s",
                    exc,
                )

    def score_prediction(self, prediction: CanonicalPrediction) -> ChemCensorScore:
        """Return official dual scores, retaining invalid and failed outputs."""
        if not prediction.is_valid:
            return ChemCensorScore(
                prediction_id=prediction.prediction.prediction_id,
                canonical_reaction_smiles=None,
                with_functional_groups=None,
                without_functional_groups=None,
                status="invalid",
                provenance=None,
                error=prediction.error or "invalid prediction",
            )

        reaction_smiles = prediction.canonical_reaction_smiles
        assert reaction_smiles is not None
        if self._cache is not None:
            cached_score = self._get_cached_score(reaction_smiles, prediction)
            if cached_score is not None:
                return cached_score
        try:
            result = self._engine.evaluate(reaction_smiles)
            with_functional_groups = _coerce_official_score(
                getattr(result, "with_functional_groups"),
                "with_functional_groups",
            )
            without_functional_groups = _coerce_official_score(
                getattr(result, "without_functional_groups"),
                "without_functional_groups",
            )
            status = _status_from_official_scores(
                with_functional_groups,
                without_functional_groups,
            )
            score = ChemCensorScore(
                prediction_id=prediction.prediction.prediction_id,
                canonical_reaction_smiles=reaction_smiles,
                with_functional_groups=with_functional_groups,
                without_functional_groups=without_functional_groups,
                status=status,
                provenance=self.provenance,
            )
        except Exception as exc:  # a per-reaction external-engine boundary
            return ChemCensorScore(
                prediction_id=prediction.prediction.prediction_id,
                canonical_reaction_smiles=reaction_smiles,
                with_functional_groups=-1.0,
                without_functional_groups=-1.0,
                status="processing_failed",
                provenance=self.provenance,
                error=f"{type(exc).__name__}: {exc}",
            )
        self._store_cached_score(score)
        return score

    def _get_cached_score(
        self,
        reaction_smiles: str,
        prediction: CanonicalPrediction,
    ) -> ChemCensorScore | None:
        """Read an optional cache without letting it alter scorer semantics."""
        assert self._cache is not None
        try:
            return self._cache.get(
                reaction_smiles,
                self.provenance,
                prediction_id=prediction.prediction.prediction_id,
            )
        except Exception as exc:  # cache is a performance optimization only
            logger.warning("ChemCensor score cache read failed; bypassing cache: %s", exc)
            return None

    def _store_cached_score(self, score: ChemCensorScore) -> None:
        """Attempt a cache write without converting a good result into failure."""
        if self._cache is None:
            return
        try:
            self._cache.put(score)
        except Exception as exc:  # cache is a performance optimization only
            logger.warning("ChemCensor score cache write failed; bypassing cache: %s", exc)


def _status_from_official_scores(
    with_functional_groups: float,
    without_functional_groups: float,
) -> str:
    """Map the official paired score sentinels to the public status vocabulary."""
    if with_functional_groups == -1.0 and without_functional_groups == -1.0:
        return "processing_failed"
    if with_functional_groups == 0.0 and without_functional_groups == 0.0:
        return "no_precedent"
    if with_functional_groups > 0.0 and without_functional_groups > 0.0:
        return "scored"
    raise ValueError("official ChemCensor result has inconsistent score variants")


def _coerce_official_score(value: Any, field_name: str) -> float:
    """Reject undocumented score shapes instead of silently changing metrics."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"official {field_name} score is not numeric")
    numeric_value = float(value)
    if not -1.0 <= numeric_value <= 5.0:
        raise ValueError(f"official {field_name} score is outside [-1, 5]")
    return numeric_value
