"""SMILES canonicalization and forward-reaction construction."""

from __future__ import annotations

from deepretro.evaluation.types import CanonicalPrediction, SingleStepPrediction
from deepretro.utils.utils_molecule import canonicalize, is_valid_smiles


def canonicalize_prediction(prediction: SingleStepPrediction) -> CanonicalPrediction:
    """Canonicalize one precursor set and construct ``precursors>>product``.

    The operation preserves duplicate precursor components because their
    multiplicity may carry experimental meaning.  It sorts only their
    canonicalized order so equivalent unordered precursor sets share a stable
    diversity key.
    """
    canonical_product = canonicalize_smiles(prediction.product_smiles)
    if canonical_product is None:
        return _invalid_prediction(
            prediction,
            f"invalid product SMILES: {prediction.product_smiles!r}",
        )

    canonical_precursors: list[str] = []
    for index, precursor in enumerate(prediction.precursor_smiles):
        canonical_precursor = canonicalize_smiles(precursor)
        if canonical_precursor is None:
            return _invalid_prediction(
                prediction,
                f"invalid precursor SMILES at index {index}: {precursor!r}",
            )
        canonical_precursors.append(canonical_precursor)

    sorted_precursors = tuple(sorted(canonical_precursors))
    reaction_smiles = f"{'.'.join(sorted_precursors)}>>{canonical_product}"
    return CanonicalPrediction(
        prediction=prediction,
        canonical_product_smiles=canonical_product,
        canonical_precursor_smiles=sorted_precursors,
        canonical_reaction_smiles=reaction_smiles,
    )


def canonicalize_smiles(smiles: str) -> str | None:
    """Return an isomeric canonical SMILES or ``None`` when invalid.

    Delegates to the package-standard helpers in
    :mod:`deepretro.utils.utils_molecule`.  ``canonicalize`` returns the input
    unchanged on a parse failure, so validity is checked first to preserve this
    function's ``None`` contract.
    """
    if not is_valid_smiles(smiles):
        return None
    return canonicalize(smiles)


def _invalid_prediction(
    prediction: SingleStepPrediction,
    error: str,
) -> CanonicalPrediction:
    """Build the common audit representation for invalid model output."""
    return CanonicalPrediction(
        prediction=prediction,
        canonical_product_smiles=None,
        canonical_precursor_smiles=None,
        canonical_reaction_smiles=None,
        error=error,
    )
