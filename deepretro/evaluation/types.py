"""Public value types for single-step plausibility evaluation."""

from __future__ import annotations

import math
import string
from dataclasses import dataclass, field
from typing import Any, Literal, Mapping


def _require_text(value: str, field_name: str) -> None:
    """Reject blank identifiers and SMILES fields at the public boundary."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")


@dataclass(frozen=True)
class SingleStepPrediction:
    """One ranked DeepRetro-style precursor-set prediction.

    The prediction represents a retrosynthetic proposal: ``product_smiles`` is
    the target-side molecule and ``precursor_smiles`` is the proposed set of
    immediate precursors.  Evaluation canonicalizes it into the forward
    reaction form required by precedent scorers.
    """

    prediction_id: str
    target_id: str
    rank: int
    product_smiles: str
    precursor_smiles: tuple[str, ...]
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate structural fields without requiring RDKit parsing yet."""
        _require_text(self.prediction_id, "prediction_id")
        _require_text(self.target_id, "target_id")
        _require_text(self.product_smiles, "product_smiles")
        if isinstance(self.rank, bool) or not isinstance(self.rank, int):
            raise ValueError("rank must be a positive integer")
        if self.rank < 1:
            raise ValueError("rank must be a positive integer")
        if not isinstance(self.precursor_smiles, tuple):
            raise ValueError("precursor_smiles must be a tuple")
        if not self.precursor_smiles:
            raise ValueError("precursor_smiles must contain at least one SMILES string")
        for index, smiles in enumerate(self.precursor_smiles):
            _require_text(smiles, f"precursor_smiles[{index}]")
        if not isinstance(self.metadata, Mapping):
            raise ValueError("metadata must be a mapping")


@dataclass(frozen=True)
class EvaluationTarget:
    """One target in the explicit benchmark denominator.

    ``reference_precursor_smiles`` is optional because exact-match remains a
    companion diagnostic rather than a requirement for plausibility scoring.
    """

    target_id: str
    product_smiles: str
    reference_precursor_smiles: tuple[str, ...] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate benchmark identity and optional single reference shape."""
        _require_text(self.target_id, "target_id")
        _require_text(self.product_smiles, "product_smiles")
        if self.reference_precursor_smiles is not None:
            if not isinstance(self.reference_precursor_smiles, tuple):
                raise ValueError("reference_precursor_smiles must be a tuple")
            if not self.reference_precursor_smiles:
                raise ValueError(
                    "reference_precursor_smiles must be non-empty when provided"
                )
            for index, smiles in enumerate(self.reference_precursor_smiles):
                _require_text(smiles, f"reference_precursor_smiles[{index}]")
        if not isinstance(self.metadata, Mapping):
            raise ValueError("metadata must be a mapping")


@dataclass(frozen=True)
class CanonicalPrediction:
    """A prediction after deterministic RDKit canonicalization.

    Invalid LLM output is represented rather than discarded.  Such a record
    has all canonical SMILES fields set to ``None`` and an explanatory error.
    """

    prediction: SingleStepPrediction
    canonical_product_smiles: str | None
    canonical_precursor_smiles: tuple[str, ...] | None
    canonical_reaction_smiles: str | None
    error: str | None = None

    def __post_init__(self) -> None:
        """Assert the valid/invalid shapes stay disjoint.

        This record is constructed only internally by
        :mod:`deepretro.evaluation.canonical`, so this is a light programmer
        guard rather than a public-input validator.
        """
        if self.canonical_reaction_smiles is None:
            assert self.error, "invalid canonical predictions require an error"
        else:
            assert (
                self.canonical_product_smiles is not None
                and self.canonical_precursor_smiles is not None
            ), "valid canonical predictions require canonical product and precursors"

    @property
    def is_valid(self) -> bool:
        """Whether this record is safe to send to an external scorer."""
        return self.canonical_reaction_smiles is not None


ScoreStatus = Literal["scored", "no_precedent", "processing_failed", "invalid"]


@dataclass(frozen=True)
class ChemCensorProvenance:
    """Versioned identity of the local official scorer/database pairing."""

    scorer_name: str
    package_version: str
    database_filename: str
    database_version: str
    database_sha256: str
    database_size_bytes: int
    max_center_type: int
    find_exact_match: bool

    def __post_init__(self) -> None:
        """Reject incomplete provenance that would make a score irreproducible."""
        for field_name in (
            "scorer_name",
            "package_version",
            "database_filename",
            "database_version",
            "database_sha256",
        ):
            _require_text(getattr(self, field_name), field_name)
        if len(self.database_sha256) != 64:
            raise ValueError("database_sha256 must be a SHA-256 hex digest")
        if any(character not in string.hexdigits for character in self.database_sha256):
            raise ValueError("database_sha256 must be a SHA-256 hex digest")
        if (
            isinstance(self.database_size_bytes, bool)
            or not isinstance(self.database_size_bytes, int)
            or self.database_size_bytes < 0
        ):
            raise ValueError("database_size_bytes must be non-negative")
        if isinstance(self.max_center_type, bool) or not isinstance(
            self.max_center_type, int
        ):
            raise ValueError("max_center_type must be an integer")
        if self.max_center_type < 1:
            raise ValueError("max_center_type must be positive")
        if not isinstance(self.find_exact_match, bool):
            raise ValueError("find_exact_match must be a boolean")


@dataclass(frozen=True)
class ChemCensorScore:
    """One normalized official ChemCensor result or a retained failure state."""

    prediction_id: str
    canonical_reaction_smiles: str | None
    with_functional_groups: float | None
    without_functional_groups: float | None
    status: ScoreStatus
    provenance: ChemCensorProvenance | None
    error: str | None = None

    def __post_init__(self) -> None:
        """Enforce status/score invariants at the public boundary."""
        _require_text(self.prediction_id, "prediction_id")
        if not isinstance(self.status, str) or self.status not in {
            "scored",
            "no_precedent",
            "processing_failed",
            "invalid",
        }:
            raise ValueError(f"unsupported ChemCensor score status: {self.status!r}")
        if self.error is not None:
            _require_text(self.error, "error")

        if self.status == "invalid":
            if self.canonical_reaction_smiles is not None:
                raise ValueError("invalid scores cannot have a canonical reaction")
            if self.provenance is not None:
                raise ValueError("invalid scores cannot have provenance")
            if (
                self.with_functional_groups is not None
                or self.without_functional_groups is not None
            ):
                raise ValueError("invalid scores cannot have numeric score values")
            return

        _require_text(self.canonical_reaction_smiles or "", "canonical_reaction_smiles")
        if self.provenance is None:
            raise ValueError("non-invalid scores require provenance")
        if not isinstance(self.provenance, ChemCensorProvenance):
            raise ValueError("provenance must be a ChemCensorProvenance instance")
        with_functional_groups = _validate_score_value(
            self.with_functional_groups,
            "with_functional_groups",
        )
        without_functional_groups = _validate_score_value(
            self.without_functional_groups,
            "without_functional_groups",
        )
        if self.status == "no_precedent" and (
            with_functional_groups != 0.0 or without_functional_groups != 0.0
        ):
            raise ValueError("no_precedent scores must be exactly 0")
        if self.status == "processing_failed" and (
            with_functional_groups != -1.0 or without_functional_groups != -1.0
        ):
            raise ValueError("processing_failed scores must be exactly -1")
        if self.status == "scored" and (
            with_functional_groups <= 0.0 or without_functional_groups <= 0.0
        ):
            raise ValueError("scored ChemCensor values must be positive")


def _validate_score_value(value: float | None, field_name: str) -> float:
    """Validate the documented ChemCensor numeric range."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field_name} must be a numeric ChemCensor score")
    numeric_value = float(value)
    if not math.isfinite(numeric_value) or not -1.0 <= numeric_value <= 5.0:
        raise ValueError(f"{field_name} must be a finite score in [-1, 5]")
    return numeric_value
