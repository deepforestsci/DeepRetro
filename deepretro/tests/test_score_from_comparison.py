"""Unit tests for the pure-arithmetic half of the hallucination checker."""

from typing import Any

import pytest

from deepretro.algorithms.hallucination_checker import (
    hallucination_compare_molecules,
    score_from_comparison,
)
from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS


def _cmp(**over: Any) -> Any:
    """A comparison dict with nothing wrong, overridden per test."""
    base = {
        "valid_reactant": True,
        "valid_product": True,
        "atom_count_consistent": True,
        "ring_size_changes": [],
        "substituent_position_changes": [],
        "detected_issues": [],
        "atom_count_deltas": {},
        "aromatic_atom_delta": 0,
        "bond_delta": 0,
        "n_ring_changes": 0,
        "n_position_changes": 0,
    }
    base.update(over)
    return base


def test_clean_scores_100_low() -> None:
    res = score_from_comparison(_cmp())
    assert res["score"] == 100
    assert res["severity"] == "low"


def test_invalid_smiles_returns_the_same_keys_as_a_normal_score() -> None:
    """Every early return carries `penalties`, so callers need no special case.

    It used to be omitted only on the invalid path, so code reading
    ``report["penalties"]`` raised KeyError on exactly the inputs most likely to
    be malformed.
    """
    res = score_from_comparison(_cmp(valid_reactant=False))
    assert sorted(res.keys()) == [
        "message",
        "penalties",
        "penalty_total",
        "score",
        "severity",
        "unassessable",
    ]
    assert res["score"] == 0
    assert res["severity"] == "critical"
    assert res["penalties"] == []


def test_atom_penalty_sums_over_elements() -> None:
    """Summed, not max: the score must match the explanation it prints.

    One "-N points" line is emitted per element, so taking max() described
    penalties that were never applied. Deltas 1 and 3 -> penalties 5 and 15 ->
    20, giving 80. Under the old max() this was 85.
    """
    assert (
        score_from_comparison(
            _cmp(atom_count_consistent=False, atom_count_deltas={"C": 1, "N": 3})
        )["score"]
        == 80
    )


def test_atom_penalty_sum_is_capped() -> None:
    """Each element is capped, and so is the total."""
    w = DEFAULT_WEIGHTS.replace(cap_atom=12)
    res = score_from_comparison(
        _cmp(atom_count_consistent=False, atom_count_deltas={"C": 1, "N": 3, "O": 2}),
        w,
    )
    assert res["score"] == 88  # total penalty clamped to cap_atom=12


def test_multi_element_case_is_actually_discriminating() -> None:
    """Guard against the golden lock's blind spot.

    The one golden fixture with several mismatched elements scores 0 under BOTH
    max() and sum(), because other penalties already floor it. So the goldens
    passed this change by accident and would not catch a regression here.
    """
    deltas = {"C": 1, "N": 3}
    summed = score_from_comparison(
        _cmp(atom_count_consistent=False, atom_count_deltas=deltas)
    )["score"]
    assert summed == 80
    assert summed != 85, "must differ from the old max() behaviour"


def test_atom_penalty_emits_one_description_per_element() -> None:
    res = score_from_comparison(
        _cmp(atom_count_consistent=False, atom_count_deltas={"C": 1, "N": 3})
    )
    assert len(res["penalties"]) == 2


def test_atom_cap() -> None:
    w = DEFAULT_WEIGHTS.replace(cap_atom=30)
    assert (
        score_from_comparison(
            _cmp(atom_count_consistent=False, atom_count_deltas={"C": 100}), w
        )["score"]
        == 70
    )


def test_ring_penalty_and_cap() -> None:
    assert score_from_comparison(_cmp(n_ring_changes=1))["score"] == 75
    assert score_from_comparison(_cmp(n_ring_changes=4))["score"] == 50


def test_position_penalty() -> None:
    assert score_from_comparison(_cmp(n_position_changes=1))["score"] == 40


def test_position_penalty_caps_and_emits_one_description() -> None:
    """Two changes hit the cap of 100, and the term emits a single line."""
    res = score_from_comparison(_cmp(n_position_changes=2))
    assert res["score"] == 0
    assert len(res["penalties"]) == 1


def test_aromaticity_fires_strictly_above_trigger() -> None:
    assert score_from_comparison(_cmp(aromatic_atom_delta=3))["score"] == 60
    assert score_from_comparison(_cmp(aromatic_atom_delta=2))["score"] == 100


def test_aromaticity_trigger_tunable_without_recomparing() -> None:
    w = DEFAULT_WEIGHTS.replace(arom_trigger=5)
    assert score_from_comparison(_cmp(aromatic_atom_delta=3), w)["score"] == 100


def test_bond_penalty_and_cap() -> None:
    assert score_from_comparison(_cmp(bond_delta=2))["score"] == 90
    assert score_from_comparison(_cmp(bond_delta=50))["score"] == 70


def test_zero_weight_disables_a_term() -> None:
    w = DEFAULT_WEIGHTS.replace(w_ring=0)
    assert score_from_comparison(_cmp(n_ring_changes=2), w)["score"] == 100


def test_score_clamped_at_zero() -> None:
    assert (
        score_from_comparison(
            _cmp(
                atom_count_consistent=False,
                atom_count_deltas={"C": 20},
                n_ring_changes=2,
                n_position_changes=2,
            )
        )["score"]
        == 0
    )


@pytest.mark.parametrize(
    "score,expected",
    [
        (100, "low"),
        (80, "low"),
        (79, "medium"),
        (40, "medium"),
        (39, "high"),
        (20, "high"),
        (19, "critical"),
        (0, "critical"),
    ],
)
def test_severity_cutoffs(score: Any, expected: Any) -> None:
    """Every boundary is reachable with w_atom=1; no case may be skipped."""
    w = DEFAULT_WEIGHTS.replace(w_atom=1, cap_atom=100)
    cmp_ = _cmp(atom_count_consistent=False, atom_count_deltas={"C": 100 - score})
    res = score_from_comparison(cmp_, w)
    assert res["score"] == score
    assert res["severity"] == expected


class _NoIssuesAccess(dict):
    """A comparison mapping that raises if ``detected_issues`` is read."""

    def __getitem__(self, key: Any) -> Any:
        if key == "detected_issues":
            raise AssertionError("score_from_comparison must not read detected_issues")
        return super().__getitem__(key)

    def get(self, key: Any, default: Any = None) -> Any:
        if key == "detected_issues":
            raise AssertionError("score_from_comparison must not read detected_issues")
        return super().get(key, default)


def test_score_half_never_touches_detected_issues() -> None:
    """The cacheability invariant, enforced by access rather than by value."""
    real = hallucination_compare_molecules("CCOCCN", "ClC(=O)c1ccccc1Br")
    guarded = _NoIssuesAccess(real)
    result = score_from_comparison(guarded)  # must not raise
    assert result["score"] == score_from_comparison(dict(real))["score"]


def test_compare_half_emits_numeric_fields() -> None:
    cmp_ = hallucination_compare_molecules("CCOCCN", "ClC(=O)c1ccccc1Br")
    assert all(
        isinstance(v, int) and v >= 0 for v in cmp_["atom_count_deltas"].values()
    )
    assert isinstance(cmp_["aromatic_atom_delta"], int)
    assert isinstance(cmp_["bond_delta"], int)
    assert cmp_["n_ring_changes"] == len(cmp_["ring_size_changes"])
    assert cmp_["n_position_changes"] == len(cmp_["substituent_position_changes"])


def test_arom_trigger_is_honoured_by_the_compare_half_too() -> None:
    """The compare half emits the issue string under the trigger it is given."""
    strict = hallucination_compare_molecules("c1ccccc1", "C1CCCCC1", arom_trigger=0)
    loose = hallucination_compare_molecules("c1ccccc1", "C1CCCCC1", arom_trigger=8)
    assert any("aromaticity" in i for i in strict["detected_issues"])
    assert not any("aromaticity" in i for i in loose["detected_issues"])
    # the numeric field is identical either way -- that is what makes the
    # comparison cacheable across trials
    assert strict["aromatic_atom_delta"] == loose["aromatic_atom_delta"]
