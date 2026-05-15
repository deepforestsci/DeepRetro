"""Unit tests for route parsing metric helpers."""

from __future__ import annotations

import inspect

from deepretro.utils.parse_metrics import ReactionMetricCalculator
from deepretro.utils.variables import REACTION_ENCODING_NAMES


class FakeReactionClassifier:
    """Minimal classifier stub returning one reaction encoding."""

    def predict(self, fingerprints: list[list[int]]) -> list[int]:
        """Return a single deterministic reaction encoding."""
        assert fingerprints == [[1, 0, 1, 0]]
        return [0]


def test_reaction_metric_calculator_returns_na_without_model_path() -> None:
    """Scalability should be unavailable when no classifier path is configured."""
    calculator = ReactionMetricCalculator(model_path="")

    assert calculator.scalability_index("CC", "CCO") == "N/A"


def test_reaction_metric_calculator_predicts_reaction_type() -> None:
    """Calculator should load the model and map classifier output to names."""
    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=lambda path: FakeReactionClassifier(),
        fingerprint_calculator=lambda smiles: [1, 0],
    )

    reaction_name, reaction_index = calculator.reaction_type("CC", "CCO")

    assert reaction_name == REACTION_ENCODING_NAMES[0]
    assert reaction_index == 0

def test_reaction_type_returns_unknown_when_fingerprint_is_none() -> None:
    """Should return UNKNOWN when fingerprint calculation fails."""

    calculator = ReactionMetricCalculator(
        model_path="model.joblib",
        model_loader=lambda path: FakeReactionClassifier(),
        fingerprint_calculator=lambda smiles: None,  # force failure
    )

    name, index = calculator.reaction_type("CC", "CCO")

    assert name == "Unknown Reaction"
    assert index == -1

def test_reaction_type_handles_missing_model_file() -> None:
    """Missing model file should return UNKNOWN_REACTION."""

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


def test_reaction_type_logs_missing_model_file() -> None:
    """Missing model files should be logged before returning UNKNOWN_REACTION."""

    class RecordingLogger:
        def __init__(self) -> None:
            self.calls = []

        def error(self, event: str, **kwargs) -> None:
            self.calls.append((event, kwargs))

    def failing_loader(path: str):
        raise FileNotFoundError(path)

    logger = RecordingLogger()
    calculator = ReactionMetricCalculator(
        model_path="missing.joblib",
        model_loader=failing_loader,
        fingerprint_calculator=lambda smiles: [1, 0],
        logger=logger,
    )

    name, index = calculator.reaction_type("CC", "CCO")

    assert name == "Unknown Reaction"
    assert index == -1
    assert logger.calls == [
        (
            "Error in metric calculation",
            {
                "function": "get_reaction_type",
                "error": "missing.joblib",
            },
        )
    ]

def test_scalability_index_returns_na_when_reaction_unknown() -> None:
    """Unknown reaction should lead to N/A scalability."""

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

