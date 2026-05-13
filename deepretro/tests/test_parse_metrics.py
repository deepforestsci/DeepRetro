"""Unit tests for route parsing metric helpers."""

from __future__ import annotations

from deepretro.utils.parse_metrics import (
    ReactionMetricCalculator,
    calc_scalability_index,
    get_reaction_type,
)
from deepretro.utils.variables import ENCODING_SCALABILITY, REACTION_ENCODING_NAMES


class DummyReactionClassifier:
    """Minimal classifier stub returning one reaction encoding."""

    def predict(self, fingerprints: list[list[int]]) -> list[int]:
        """Return a single deterministic reaction encoding."""
        assert fingerprints == [[1, 0, 1, 0]]
        return [0]


def test_reaction_metric_calculator_returns_na_without_model_path() -> None:
    """When model_path is empty, no classifier is loaded and scalability_index
    should return 'N/A' for any molecule pair."""
    calculator = ReactionMetricCalculator(model_path="")

    assert calculator.scalability_index("CC", "CCO") == "N/A"


def test_reaction_metric_calculator_predicts_reaction_type() -> None:
    """Happy path: with an injected model and fingerprint calculator, the
    calculator should load the classifier, concatenate fingerprints [1,0]+[1,0],
    and map the prediction index (0) to the corresponding reaction name."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=lambda path: DummyReactionClassifier(),
        fingerprint_calculator=lambda smiles: [1, 0],
    )

    reaction_name, reaction_index = calculator.reaction_type("CC", "CCO")

    assert reaction_name == REACTION_ENCODING_NAMES[0]
    assert reaction_index == 0


def test_reaction_type_returns_unknown_when_fingerprint_is_none() -> None:
    """When the fingerprint calculator returns None (e.g. invalid SMILES or
    RDKit failure), reaction_type should gracefully return UNKNOWN_REACTION
    instead of crashing."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=lambda path: DummyReactionClassifier(),
        fingerprint_calculator=lambda smiles: None,
    )

    name, index = calculator.reaction_type("CC", "CCO")

    assert name == "Unknown Reaction"
    assert index == -1


def test_reaction_type_handles_missing_model_file() -> None:
    """When the model loader raises FileNotFoundError (model file missing
    from disk), reaction_type should catch it and return UNKNOWN_REACTION
    so that route parsing can continue without a classifier."""

    def failing_loader(path: str):
        raise FileNotFoundError

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=failing_loader,
        fingerprint_calculator=lambda smiles: [1, 0],
    )

    name, index = calculator.reaction_type("CC", "CCO")

    assert name == "Unknown Reaction"
    assert index == -1


def test_scalability_index_returns_na_when_reaction_unknown() -> None:
    """When the classifier returns an index (999) that doesn't map to a known
    reaction encoding, scalability_index should return 'N/A' rather than
    raising a KeyError."""

    class BadClassifier:
        def predict(self, x):
            return [999]

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=lambda path: BadClassifier(),
        fingerprint_calculator=lambda smiles: [1, 0],
    )

    result = calculator.scalability_index("CC", "CCO")

    assert result == "N/A"


def test_scalability_index_returns_label_for_valid_prediction() -> None:
    """Happy path: when the classifier predicts index 0, scalability_index
    should look up and return the corresponding label from ENCODING_SCALABILITY,
    verifying the full prediction-to-label pipeline."""

    class PredictZero:
        def predict(self, x):
            return [0]

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=lambda path: PredictZero(),
        fingerprint_calculator=lambda smiles: [1, 0],
    )

    result = calculator.scalability_index("CC", "CCO")

    assert result == str(ENCODING_SCALABILITY[0])


def test_reaction_type_logs_error_on_file_not_found() -> None:
    """FileNotFoundError should be logged via the injected logger before
    returning UNKNOWN_REACTION. Verifies the fix for silently swallowed
    errors that made missing model files hard to debug."""
    logged = []

    class CapturingLogger:
        def error(self, *args, **kwargs):
            logged.append((args, kwargs))

    def failing_loader(path: str):
        raise FileNotFoundError("model.joblib not found")

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=failing_loader,
        fingerprint_calculator=lambda smiles: [1, 0],
        logger=CapturingLogger(),
    )

    calculator.reaction_type("CC", "CCO")

    assert len(logged) == 1
    assert "model.joblib not found" in logged[0][1]["error"]


def test_get_reaction_type_compatibility_wrapper() -> None:
    """The module-level get_reaction_type() backward-compatible wrapper should
    delegate to ReactionMetricCalculator and return the same result."""
    result = get_reaction_type("CC", "CCO", model_path="")

    assert result == ("Unknown Reaction", -1)


def test_calc_scalability_index_compatibility_wrapper() -> None:
    """The module-level calc_scalability_index() backward-compatible wrapper
    should delegate to ReactionMetricCalculator using the env-configured
    model path (empty by default in tests), returning 'N/A'."""
    result = calc_scalability_index("CC", "CCO")

    assert result == "N/A"
