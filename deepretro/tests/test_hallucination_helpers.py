"""Tests for hallucination checker adapter helpers."""

from __future__ import annotations

from unittest.mock import MagicMock
from unittest.mock import call

import pytest

from deepretro.models.hallucination_helpers import (
    build_ml_checker,
    filter_with_checker,
    resolve_hallucination,
)


def test_build_ml_checker_keeps_non_hallucinated_pathways(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """ML checker keeps pathways whose classifier prediction is clean."""
    classifier = MagicMock()
    classifier.predict_single.side_effect = [
        {"is_hallucination": True},
        {"is_hallucination": False},
    ]
    monkeypatch.setattr(
        "deepretro.utils.utils_molecule.is_valid_smiles",
        lambda smiles: True,
    )

    status, retained = build_ml_checker(classifier)("CCO", [["CC"], ["CO"]])

    assert status == 200
    assert retained == [["CO"]]
    classifier.predict_single.assert_has_calls(
        [
            call("CCO", "CC"),
            call("CCO", "CO"),
        ]
    )


def test_build_ml_checker_drops_invalid_reactant_sets(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Invalid joined reactant SMILES are ignored before classifier calls."""
    classifier = MagicMock()
    monkeypatch.setattr(
        "deepretro.utils.utils_molecule.is_valid_smiles",
        lambda smiles: False,
    )

    status, retained = build_ml_checker(classifier)("CCO", [["not-valid"]])

    assert status == 200
    assert retained == []
    classifier.predict_single.assert_not_called()


def test_filter_with_checker_preserves_explanation_and_confidence_alignment() -> None:
    """Filtering retains metadata aligned to the surviving pathway."""

    def checker(product: str, pathways: list[list[str]]) -> tuple[int, list[list[str]]]:
        assert product == "CCO"
        assert pathways == [["CC"], ["CO"]]
        return 200, [["CO"]]

    pathways, explanations, confidence = filter_with_checker(
        "CCO",
        [["CC"], ["CO"]],
        ["first", "second"],
        [0.1, 0.9],
        checker,
    )

    assert pathways == [["CO"]]
    assert explanations == ["second"]
    assert confidence == [0.9]


def test_filter_with_checker_short_circuits_without_checker() -> None:
    """No checker leaves pathways and metadata unchanged."""
    pathways, explanations, confidence = filter_with_checker(
        "CCO",
        [["CC"], ["CO"]],
        ["first", "second"],
        [0.1, 0.9],
        None,
    )

    assert pathways == [["CC"], ["CO"]]
    assert explanations == ["first", "second"]
    assert confidence == [0.1, 0.9]


def test_filter_with_checker_returns_empty_on_checker_failure() -> None:
    """Non-200 checker status drops all aligned metadata."""

    def checker(product: str, pathways: list[list[str]]) -> tuple[int, list[list[str]]]:
        return 500, pathways

    assert filter_with_checker("CCO", [["CC"]], ["first"], [0.1], checker) == (
        [],
        [],
        [],
    )


def test_filter_with_checker_handles_duplicate_pathways_in_order() -> None:
    """Duplicate retained pathways consume metadata entries one at a time."""

    def checker(product: str, pathways: list[list[str]]) -> tuple[int, list[list[str]]]:
        return 200, [["CC"], ["CC"]]

    pathways, explanations, confidence = filter_with_checker(
        "CCO",
        [["CC"], ["CC"], ["CO"]],
        ["first", "second", "third"],
        [0.1, 0.2, 0.3],
        checker,
    )

    assert pathways == [["CC"], ["CC"]]
    assert explanations == ["first", "second"]
    assert confidence == [0.1, 0.2]


def test_resolve_hallucination_modes() -> None:
    """Mode resolution is explicit and validates ML classifier inputs."""
    assert resolve_hallucination("none", None) is None
    assert resolve_hallucination("heuristic", None) is None

    classifier = MagicMock()
    classifier.predict_single.return_value = {"is_hallucination": False}
    assert callable(resolve_hallucination("ml", classifier))

    with pytest.raises(ValueError, match="hallucination_mode"):
        resolve_hallucination("unknown", None)
    with pytest.raises(ValueError, match="requires"):
        resolve_hallucination("ml", None)
