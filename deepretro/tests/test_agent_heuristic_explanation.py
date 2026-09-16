"""Regression tests for heuristic feedback delivered through the agent registry."""

from __future__ import annotations

import json

import pytest

from deepretro.agents.tools import build_tool_registry
from deepretro.algorithms.hallucination_checker import hallucination_compare_molecules
from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS
from deepretro.models.hallucination_checker import HallucinationChecker
from deepretro.models.hallucination_helpers import resolve_hallucination


@pytest.mark.parametrize("resolved", [False, True])
@pytest.mark.parametrize("reactants", ["CC", ["CC"]])
def test_retained_flagged_step_reaches_agent_with_explanation(
    resolved: bool, reactants: str | list[str]
) -> None:
    """Both entry points preserve details and flag a retained fallback."""
    checker = (
        resolve_hallucination("heuristic", None)
        if resolved
        else HallucinationChecker(checker_type="heuristic")
    )
    assert checker is not None
    assert checker("c1ccccc1", [["CC"]]) == (200, [["CC"]])
    registry = build_tool_registry(hallucination_checker=checker)
    result = registry.execute(
        "check_hallucination", {"product": "c1ccccc1", "reactants": reactants}
    )
    assert result["is_hallucination"] is True
    assert result["score"] == 0
    assert result["severity"] == "critical"
    assert result["unassessable"] is False
    assert "detected_issues" not in result
    assert (
        "Atom count mismatch for C: Reactant has 2, Product has 6"
        in result["explanation"]["detected_issues"]
    )
    assert "6-membered ring added" in result["explanation"]["ring_size_changes"]
    comparison = hallucination_compare_molecules("CC", "c1ccccc1")
    for key, value in result["explanation"].items():
        assert value == comparison[key]
    assert json.loads(json.dumps(result)) == result


def test_clean_step_has_empty_issues() -> None:
    """An unflagged, assessable transformation retains the explanation shape."""
    checker = HallucinationChecker(checker_type="heuristic")
    result = build_tool_registry(hallucination_checker=checker).execute(
        "check_hallucination", {"product": "CC=O", "reactants": ["CCO"]}
    )
    assert result["is_hallucination"] is False
    assert result["score"] == 100
    assert result["explanation"]["detected_issues"] == []


@pytest.mark.parametrize(
    ("product", "reactants", "reason"),
    [
        ("CCO", "invalid", "Invalid"),
        ("invalid", "CCO", "Invalid"),
        ("CCO", "", "Invalid"),
        ("CCO", [], "Invalid"),
        ("CCO", "CCO", "Degenerate"),
        ("CCO", "[CH2]CO", "radical"),
    ],
)
def test_unassessable_step_has_nullable_verdict(
    product: str, reactants: str | list[str], reason: str
) -> None:
    """Unsupported inputs convey their reason without a chemical verdict."""
    checker = HallucinationChecker(checker_type="heuristic")
    result = build_tool_registry(hallucination_checker=checker).execute(
        "check_hallucination", {"product": product, "reactants": reactants}
    )
    assert result["is_hallucination"] is None
    assert result["unassessable"] is True
    assert reason in result["message"]
    assert isinstance(result["explanation"]["detected_issues"], list)


def test_assessment_uses_custom_threshold_and_weights() -> None:
    """Severity and retention cannot override the configured flag threshold."""
    weights = DEFAULT_WEIGHTS.replace(w_arom=0, reject_below=0)
    checker = HallucinationChecker(checker_type="heuristic", weights=weights)
    result = checker.assess("c1ccccc1", "CC")
    assert result["score"] == 30
    assert result["flagged"] is False
    assert checker.check_single_pathway("c1ccccc1", "CC") == 0


def test_substituent_details_are_preserved() -> None:
    """Relative substituent positions reach the agent without being summarized."""
    reactants, product = "Cc1ccc(O)cc1", "Cc1cccc(O)c1"
    checker = HallucinationChecker(checker_type="heuristic")
    result = build_tool_registry(hallucination_checker=checker).execute(
        "check_hallucination", {"product": product, "reactants": [reactants]}
    )
    expected = hallucination_compare_molecules(reactants, product)
    assert expected["substituent_position_changes"]
    assert (
        result["explanation"]["substituent_position_changes"]
        == expected["substituent_position_changes"]
    )


@pytest.mark.parametrize("prediction", [0, 1])
def test_ml_assessment_preserves_prediction(
    monkeypatch: pytest.MonkeyPatch, prediction: int
) -> None:
    """The new assessment path leaves ML decisions to the existing classifier."""
    checker = HallucinationChecker(checker_type="heuristic")
    checker.checker_type = "ml"
    monkeypatch.setattr(checker, "_check_pathway_ml", lambda *_: prediction)
    result = build_tool_registry(hallucination_checker=checker).execute(
        "check_hallucination", {"product": "CC=O", "reactants": ["CCO"]}
    )
    assert result["is_hallucination"] is bool(prediction)
    assert result["source"] == "ml"
    assert "explanation" not in result
