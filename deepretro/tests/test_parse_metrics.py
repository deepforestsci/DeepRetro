"""Unit tests for route parsing metric helpers."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, NoReturn, Optional

import structlog.testing

from deepretro.utils.parse_metrics import (
    ReactionMetricCalculator,
    calc_scalability_index,
    get_reaction_type,
)
from deepretro.utils.variables import ENCODING_SCALABILITY, REACTION_ENCODING_NAMES


def test_reaction_metric_calculator_returns_na_without_model_path() -> None:
    """When model_path is empty, no classifier is loaded and scalability_index
    should return 0 for any molecule pair."""
    calculator = ReactionMetricCalculator(model_path="")

    assert calculator.scalability_index("CC", "CCO") == 0


def test_reaction_metric_calculator_predicts_reaction_type(
    dummy_model_loader: Callable[[str], Any],
    fingerprint_calculator: Callable[[str], list[int]],
) -> None:
    """Happy path: with an injected model and fingerprint calculator, the
    calculator should load the classifier, concatenate fingerprints [1,0]+[1,0],
    and map the prediction index (0) to the corresponding reaction name."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=dummy_model_loader,
        fingerprint_calculator=fingerprint_calculator,
    )

    reaction_name, reaction_index = calculator.reaction_type("CC", "CCO")

    assert reaction_name == REACTION_ENCODING_NAMES[0]
    assert reaction_index == 0


def test_reaction_type_returns_unknown_when_fingerprint_is_none(
    dummy_model_loader: Callable[[str], Any],
    missing_fingerprint_calculator: Callable[[str], Optional[list[int]]],
) -> None:
    """When the fingerprint calculator returns None (e.g. invalid SMILES or
    RDKit failure), reaction_type should gracefully return UNKNOWN_REACTION
    instead of crashing."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=dummy_model_loader,
        fingerprint_calculator=missing_fingerprint_calculator,
    )

    name, index = calculator.reaction_type("CC", "CCO")

    assert name == "Unknown Reaction"
    assert index == -1


def test_reaction_type_handles_missing_model_file(
    missing_model_loader: Callable[[str], NoReturn],
    fingerprint_calculator: Callable[[str], list[int]],
) -> None:
    """When the model loader raises FileNotFoundError (model file missing
    from disk), reaction_type should catch it and return UNKNOWN_REACTION
    so that route parsing can continue without a classifier."""

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=missing_model_loader,
        fingerprint_calculator=fingerprint_calculator,
    )

    name, index = calculator.reaction_type("CC", "CCO")

    assert name == "Unknown Reaction"
    assert index == -1


def test_scalability_index_returns_zero_when_reaction_unknown(
    unknown_model_loader: Callable[[str], Any],
    fingerprint_calculator: Callable[[str], list[int]],
) -> None:
    """When the classifier returns an index (999) that doesn't map to a known
    reaction encoding, scalability_index should return 0 rather than
    raising a KeyError."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=unknown_model_loader,
        fingerprint_calculator=fingerprint_calculator,
    )

    result = calculator.scalability_index("CC", "CCO")

    assert result == 0


def test_scalability_index_returns_score_for_valid_prediction(
    dummy_model_loader: Callable[[str], Any],
    fingerprint_calculator: Callable[[str], list[int]],
) -> None:
    """Happy path: when the classifier predicts index 0, scalability_index
    should look up and return the corresponding score from ENCODING_SCALABILITY,
    verifying the full prediction-to-score pipeline."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=dummy_model_loader,
        fingerprint_calculator=fingerprint_calculator,
    )

    result = calculator.scalability_index("CC", "CCO")

    assert result == ENCODING_SCALABILITY[0]


def test_reaction_type_logs_error_on_file_not_found(
    missing_model_loader_with_message: Callable[[str], NoReturn],
    fingerprint_calculator: Callable[[str], list[int]],
) -> None:
    """FileNotFoundError should be logged via structlog before returning
    UNKNOWN_REACTION. Verifies the fix for silently swallowed errors that
    made missing model files hard to debug."""

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=missing_model_loader_with_message,
        fingerprint_calculator=fingerprint_calculator,
    )

    with structlog.testing.capture_logs() as captured:
        calculator.reaction_type("CC", "CCO")

    assert len(captured) == 1
    assert captured[0]["function"] == "get_reaction_type"
    assert "model.joblib not found" in captured[0]["error"]


def test_get_reaction_type_compatibility_wrapper() -> None:
    """The module-level get_reaction_type() backward-compatible wrapper should
    delegate to ReactionMetricCalculator and return the same result."""
    result = get_reaction_type("CC", "CCO", model_path="")

    assert result == ("Unknown Reaction", -1)


def test_calc_scalability_index_compatibility_wrapper() -> None:
    """The module-level calc_scalability_index() backward-compatible wrapper
    should delegate to ReactionMetricCalculator using the env-configured
    model path (empty by default in tests), returning 0."""
    result = calc_scalability_index("CC", "CCO")

    assert result == 0
