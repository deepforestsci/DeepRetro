"""Metric helpers used while formatting retrosynthesis route trees."""

from __future__ import annotations

import os
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, Optional

import joblib
import structlog

from deepretro.utils.variables import ENCODING_SCALABILITY, REACTION_ENCODING_NAMES

ENABLE_LOGGING = os.getenv("ENABLE_LOGGING", "true").lower() != "false"
UNKNOWN_REACTION = ("Unknown Reaction", -1)
logger = structlog.get_logger(__name__)


def _repo_root() -> Path:
    """Locate the repository root from this module path."""
    here = Path(__file__).resolve()
    for directory in [here, *here.parents]:
        if (directory / "config" / "langfuse_config.json").is_file():
            return directory
    raise FileNotFoundError(
        f"Could not locate repo root from {here} (expected config/langfuse_config.json)"
    )


def _rxn_classification_model_path() -> str:
    """Resolve the optional reaction classification model path."""
    configured_path = os.getenv("RXN_CLASSIFICATION_MODEL_PATH", "") or ""
    if not configured_path.strip():
        return ""

    path = Path(configured_path)
    return str(path if path.is_absolute() else _repo_root() / path)


def _default_fingerprint_calculator(smiles: str) -> Optional[list[int]]:
    """Calculate a model fingerprint using the molecule utility module."""
    from deepretro.utils import utils_molecule

    return utils_molecule.compute_fingerprint(smiles)


class ReactionMetricCalculator:
    """
    Calculate scalability metrics for parsed route steps.

    Examples
    --------
    >>> calculator = ReactionMetricCalculator(model_path="")
    >>> calculator.scalability_index("CC", "CCO")
    'N/A'
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        model_loader: Callable[[str], Any] = joblib.load,
        fingerprint_calculator: Callable[[str], Optional[list[int]]] = (
            _default_fingerprint_calculator
        ),
        reaction_encoding_names: Mapping[int, str] = REACTION_ENCODING_NAMES,
        scalability_encoding: Mapping[int, Any] = ENCODING_SCALABILITY,
        logger: Optional[Any] = None,
    ) -> None:
        """
        Parameters
        ----------
        model_path : str, optional
            Path to the joblib reaction classifier. ``None`` means resolve from
            ``RXN_CLASSIFICATION_MODEL_PATH``; an empty string disables classifier
            backed scalability.
        model_loader : Callable[[str], Any], optional
            Loader used to deserialize the classifier. Defaults to ``joblib.load``.
        fingerprint_calculator : Callable[[str], list[int] | None], optional
            Function used to convert SMILES to model fingerprints.
        reaction_encoding_names : Mapping[int, str], optional
            Mapping from classifier output index to reaction type label.
        scalability_encoding : Mapping[int, Any], optional
            Mapping from classifier output index to scalability label.
        logger : object, optional
            Logger exposing an ``error`` method for recoverable metric errors.
        """
        self.model_path = (
            _rxn_classification_model_path() if model_path is None else model_path
        )
        self.model_loader = model_loader
        self.fingerprint_calculator = fingerprint_calculator
        self.reaction_encoding_names = reaction_encoding_names
        self.scalability_encoding = scalability_encoding
        self.logger = logger

    def reaction_type(self, mol1: str, mol2: str) -> tuple[str, int]:
        """
        Predict the reaction type for a reactant/product pair.

        Parameters
        ----------
        mol1 : str
            Reactant SMILES.
        mol2 : str
            Product SMILES.

        Returns
        -------
        tuple[str, int]
            Reaction label and classifier encoding index. Unknown or failed
            predictions return ``("Unknown Reaction", -1)``.
        """
        try:
            return self._predict_reaction_type(mol1, mol2)
        except FileNotFoundError:
            return UNKNOWN_REACTION
        except Exception as exc:
            self._log_exception("get_reaction_type", exc)
            return UNKNOWN_REACTION

    def scalability_index(self, mol1: str, mol2: str) -> str:
        """
        Calculate the scalability label for a reactant/product pair.

        Parameters
        ----------
        mol1 : str
            Reactant SMILES.
        mol2 : str
            Product SMILES.

        Returns
        -------
        str
            Scalability label, or ``"N/A"`` when the classifier is unavailable
            or cannot classify the pair.
        """
        if not self.model_path:
            return "N/A"

        try:
            _, reaction_index = self.reaction_type(mol1, mol2)
            if reaction_index == -1:
                return "N/A"
            return str(self.scalability_encoding[reaction_index])
        except Exception as exc:
            self._log_exception("calc_scalability_index", exc)
            return "N/A"

    def _predict_reaction_type(self, mol1: str, mol2: str) -> tuple[str, int]:
        """Run the classifier for a reactant/product SMILES pair."""
        model = self.model_loader(self.model_path)
        mol1_fingerprint = self.fingerprint_calculator(mol1)
        mol2_fingerprint = self.fingerprint_calculator(mol2)
        if mol1_fingerprint is None or mol2_fingerprint is None:
            return UNKNOWN_REACTION

        prediction = model.predict([mol1_fingerprint + mol2_fingerprint])
        reaction_index = int(prediction[0])
        return self.reaction_encoding_names[reaction_index], reaction_index

    def _log_exception(self, function_name: str, exc: Exception) -> None:
        """Log a recoverable metric calculation error."""
        active_logger = self.logger
        if active_logger is None and ENABLE_LOGGING:
            active_logger = logger
        if active_logger is not None:
            active_logger.error(
                "Error in metric calculation",
                function=function_name,
                error=str(exc),
            )


def get_reaction_type(mol1: str, mol2: str, model_path: str) -> tuple[str, int]:
    """
    Predict reaction type for a reactant/product pair.

    This compatibility wrapper delegates to ``ReactionMetricCalculator``.
    """
    return ReactionMetricCalculator(model_path=model_path).reaction_type(mol1, mol2)


def calc_scalability_index(mol1: str, mol2: str) -> str:
    """
    Calculate the scalability label for a reactant/product pair.

    This compatibility wrapper resolves the configured classifier path from
    ``RXN_CLASSIFICATION_MODEL_PATH``.
    """
    return ReactionMetricCalculator().scalability_index(mol1, mol2)
