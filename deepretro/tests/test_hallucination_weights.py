"""Unit tests for the tunable hallucination-checker weights."""

import json
from dataclasses import asdict, fields
from typing import Any

import pytest

from deepretro.algorithms.hallucination_weights import (
    DEFAULT_WEIGHTS,
    HallucinationWeights,
)


def test_defaults_match_current_hardcoded_values() -> None:
    w = HallucinationWeights()
    assert (w.w_atom, w.cap_atom) == (5, 100)
    assert (w.w_ring, w.cap_ring) == (25, 50)
    assert (w.w_pos, w.cap_pos) == (60, 100)
    assert (w.w_arom, w.arom_trigger) == (40, 2)
    assert (w.w_bond, w.cap_bond) == (5, 30)
    assert (w.cut_low, w.cut_medium, w.cut_high) == (80, 40, 20)
    assert w == DEFAULT_WEIGHTS


def test_frozen() -> None:
    with pytest.raises(Exception):
        DEFAULT_WEIGHTS.w_atom = 7


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cut_low": 30, "cut_medium": 40, "cut_high": 20},  # medium >= low
        {"cut_medium": 10},  # medium <= high
        {"cut_medium": 20},  # medium == high
        {"cut_medium": 80},  # medium == low
    ],
)
def test_cutoff_ordering_enforced(kwargs: Any) -> None:
    with pytest.raises(ValueError, match="cut_high < cut_medium < cut_low"):
        HallucinationWeights(**kwargs)


def test_negative_weight_rejected() -> None:
    with pytest.raises(ValueError, match="must be non-negative"):
        HallucinationWeights(w_atom=-1)


@pytest.mark.parametrize("bad", [1.5, "5", None, True])
def test_non_int_weight_rejected(bad: Any) -> None:
    """TypeError, not ValueError: a float weight is a type error, and bool is
    excluded explicitly because bool is a subclass of int."""
    with pytest.raises(TypeError, match="must be an int"):
        HallucinationWeights(w_atom=bad)


def test_json_round_trip(tmp_path: Any) -> None:
    w = HallucinationWeights(w_atom=9, cap_ring=77, cut_medium=55)
    path = tmp_path / "w.json"
    w.to_json(path)
    assert HallucinationWeights.from_json(path) == w
    assert len(json.loads(path.read_text())) == len(fields(w))


def test_from_json_rejects_unknown_fields(tmp_path: Any) -> None:
    path = tmp_path / "w.json"
    path.write_text(json.dumps({**asdict(DEFAULT_WEIGHTS), "extra": 1}))
    with pytest.raises(ValueError, match="unknown weight fields"):
        HallucinationWeights.from_json(path)


def test_from_json_rejects_missing_fields(tmp_path: Any) -> None:
    """best_weights.json must hold all configuration fields, not a subset."""
    path = tmp_path / "w.json"
    partial = asdict(DEFAULT_WEIGHTS)
    del partial["w_atom"]
    path.write_text(json.dumps(partial))
    with pytest.raises(ValueError, match="missing weight fields"):
        HallucinationWeights.from_json(path)


def test_replace_returns_new_instance() -> None:
    w2 = DEFAULT_WEIGHTS.replace(w_atom=11)
    assert w2.w_atom == 11 and DEFAULT_WEIGHTS.w_atom == 5


def test_reject_below_defaults_to_the_historical_gate() -> None:
    """40 is exactly the old cut_medium, so behaviour is unchanged by default."""
    assert DEFAULT_WEIGHTS.reject_below == 40
    assert DEFAULT_WEIGHTS.reject_below == DEFAULT_WEIGHTS.cut_medium


def test_reject_below_is_free_of_the_cutoff_ordering_invariant() -> None:
    """The point of the field: the gate can now reach values cut_medium cannot.

    cut_high < cut_medium forced the tuner's range to start at 21, so "almost
    never reject" was unreachable.
    """
    for value in (0, 1, 20, 101):
        assert DEFAULT_WEIGHTS.replace(reject_below=value).reject_below == value
    with pytest.raises(ValueError, match="reject_below must be <= 101"):
        DEFAULT_WEIGHTS.replace(reject_below=102)


def test_old_weight_files_without_reject_below_keep_their_meaning(
    tmp_path: Any,
) -> None:
    """Backward compatibility: files saved before the field existed.

    They gated on cut_medium, so reject_below must default to that file's
    cut_medium rather than to today's default of 40.
    """
    import json

    path = tmp_path / "old.json"
    payload = {
        f.name: getattr(DEFAULT_WEIGHTS, f.name)
        for f in fields(DEFAULT_WEIGHTS)
        if f.name != "reject_below"
    }
    payload["cut_medium"] = 71
    path.write_text(json.dumps(payload))

    loaded = HallucinationWeights.from_json(path)
    assert loaded.reject_below == 71, "must inherit the file's own gate, not 40"
