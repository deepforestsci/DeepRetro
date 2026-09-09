"""Public API regressions for structural heuristic v2."""

from typing import Any

import pytest

from deepretro.algorithms.hallucination_checker import (
    calculate_hallucination_score,
    hallucination_compare_molecules,
)
from deepretro.algorithms.hallucination_weights import HallucinationWeights
from deepretro.models.hallucination_checker import HallucinationChecker


@pytest.mark.parametrize("reactants,product", [("", "CCO"), ("CCO", ""), ("", "")])
def test_empty_structures_are_unassessable(reactants: str, product: str) -> None:
    assert calculate_hallucination_score(reactants, product)["unassessable"]


def test_map_metadata_does_not_make_a_transformation() -> None:
    report = calculate_hallucination_score("[CH3:1][CH2:2][OH:3]", "CCO")
    assert report["unassessable"]
    assert report["score"] == 0


@pytest.mark.parametrize("reverse", [False, True])
def test_ring_changes_count_only_the_difference(reverse: bool) -> None:
    pair = ["C1CCCCC1.C1CCCCC1", "C1CCCCC1"]
    if reverse:
        pair.reverse()
    report = hallucination_compare_molecules(*pair, strip_reagents=False)
    assert report["n_ring_changes"] == 1


def test_fragment_permutation_preserves_ring_matching() -> None:
    a, b = "Cc1ccccc1Cl", "Cc1cccc(Cl)c1"
    reports = [
        calculate_hallucination_score(r, p)
        for r in [f"{a}.{b}", f"{b}.{a}"]
        for p in [f"{a}.{b}.O", f"{b}.{a}.O", f"O.{a}.{b}"]
    ]
    assert all(r == reports[0] for r in reports)
    assert reports[0]["score"] == 95


def test_retained_fallback_is_still_flagged() -> None:
    checker = HallucinationChecker(
        checker_type="heuristic", weights=HallucinationWeights(reject_below=101)
    )
    assert checker("CC=O", [["CCO"]])[1] == [["CCO"]]
    assert checker.check_single_pathway("CC=O", "CCO") == 1


def test_weight_file_requires_json_object(tmp_path: Any) -> None:
    path = tmp_path / "weights.json"
    path.write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="object"):
        HallucinationWeights.from_json(path)
