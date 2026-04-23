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


def test_private_metric_calculator_helpers_are_documented_and_typed() -> None:
    """Private metric helpers should still have docstrings and type annotations."""
    private_helpers = [
        member
        for name, member in inspect.getmembers(
            ReactionMetricCalculator, inspect.isfunction
        )
        if name.startswith("_") and not name.startswith("__")
    ]

    assert private_helpers
    for helper in private_helpers:
        signature = inspect.signature(helper)
        assert inspect.getdoc(helper), helper.__name__
        assert all(
            parameter.annotation is not inspect.Signature.empty
            for parameter in signature.parameters.values()
            if parameter.name not in {"self", "cls"}
        ), helper.__name__
        assert signature.return_annotation is not inspect.Signature.empty, (
            helper.__name__
        )


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
