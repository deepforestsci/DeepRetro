"""Tests for hallucination_mode='none'."""

from typing import Any

import pytest

from deepretro.algorithms.autosolve import AutoSolver
from deepretro.models.hallucination_helpers import filter_with_checker


def test_none_mode_constructs_and_disables_the_checker() -> None:
    """It is documented in the class docstring and used in 8 doctests.

    Before this fix it raised ValueError, so the gate could never be switched
    off -- meaning nobody could measure whether the checker helps at all.
    """
    solver = AutoSolver(hallucination_mode="none")
    assert solver.hallucination_checker is None
    assert solver.hallucination_mode == "none"


def test_heuristic_and_ml_still_build_a_checker() -> None:
    solver = AutoSolver(hallucination_mode="heuristic")
    assert solver.hallucination_checker is not None


def test_mode_is_case_insensitive() -> None:
    assert AutoSolver(hallucination_mode="NONE").hallucination_checker is None


def test_unknown_mode_names_the_allowed_values() -> None:
    with pytest.raises(ValueError, match="Allowed: 'heuristic', 'ml', 'none'"):
        AutoSolver(hallucination_mode="bogus")


def test_none_checker_passes_every_pathway_through() -> None:
    """The consumer contract that makes the None checker safe."""
    pathways = [["CCO"], ["CCN"]]
    kept, expl, conf = filter_with_checker(
        "CC=O", pathways, ["a", "b"], [0.9, 0.8], None
    )
    assert kept == pathways
    assert expl == ["a", "b"]
    assert conf == [0.9, 0.8]


def test_tuned_weights_still_reach_the_heuristic_checker() -> None:
    """Guard the trap: resolve_hallucination() drops weights for heuristic mode.

    The fix must not route through it, or --hallucination-weights would be
    silently ignored on every run.
    """
    from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS

    tuned = DEFAULT_WEIGHTS.replace(cut_medium=71)
    solver = AutoSolver(hallucination_mode="heuristic", hallucination_weights=tuned)
    assert solver.hallucination_checker.weights is tuned
    assert solver.hallucination_checker.weights.cut_medium == 71


# --- The verdict must survive onto the route ---------------------------------


def test_verdict_is_recorded_on_the_route_step() -> None:
    """Ranking lets a flagged step reach a finished route; it must be visible.

    Without this the chemist gets a silently degraded route: the checker
    computed a verdict and threw it away.
    """
    from deepretro.algorithms.autosolve import reaction_tree, unsolved_leaf

    solver = AutoSolver(hallucination_mode="heuristic")
    verdict = solver._verdict_for(["CCO"], "c1ccc2ccccc2c1")
    assert verdict["flagged"] is True
    assert verdict["score"] == 0

    tree = reaction_tree("c1ccc2ccccc2c1", [unsolved_leaf("CCO")], [0.9], verdict)
    stored = tree["children"][0]["metadata"]["hallucination"]
    assert stored["flagged"] is True
    assert stored["severity"] == "critical"


def test_a_clean_step_is_recorded_as_unflagged() -> None:
    solver = AutoSolver(hallucination_mode="heuristic")
    verdict = solver._verdict_for(["CCO"], "CC=O")
    assert verdict == {
        "source": "heuristic",
        "score": 100,
        "severity": "low",
        "flagged": False,
    }


def test_no_verdict_when_the_checker_is_disabled() -> None:
    assert AutoSolver(hallucination_mode="none")._verdict_for(["CCO"], "CC=O") is None


def test_reaction_tree_omits_the_key_when_there_is_no_verdict() -> None:
    """Routes from a gate-off run must not carry an empty hallucination key."""
    from deepretro.algorithms.autosolve import reaction_tree, unsolved_leaf

    tree = reaction_tree("CC=O", [unsolved_leaf("CCO")], [0.9])
    assert "hallucination" not in tree["children"][0]["metadata"]


def test_verdict_uses_the_solvers_tuned_weights() -> None:
    """A run with tuned weights must record the verdict those weights give."""
    from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS

    strict = DEFAULT_WEIGHTS.replace(reject_below=101)
    solver = AutoSolver(hallucination_mode="heuristic", hallucination_weights=strict)
    # Everything is below a threshold of 101, so even a perfect step is flagged.
    assert solver._verdict_for(["CCO"], "CC=O")["flagged"] is True


def test_verdict_never_raises_on_bad_input() -> None:
    """Annotation is best-effort; it must not break a solve."""
    solver = AutoSolver(hallucination_mode="heuristic")
    assert solver._verdict_for([], "CC=O") is None
    assert solver._verdict_for(["not_a_smiles!!"], "CC=O") is not None


# --- An unparseable target must not become a "delivered route" ---------------


def test_unparseable_target_writes_an_error_not_a_fake_route(tmp_path: Any) -> None:
    """Three of the 14 benchmark targets were in this state.

    They ran the whole pipeline and emitted a one-step pathway with an empty
    reactant list, which the report counted as a delivered route -- so
    "solved 7 of 14" was measured against a denominator containing three targets
    nothing could ever solve.
    """
    import json

    from deepretro.batch import run_batch

    calls = []

    def fake_solve(smiles: Any) -> Any:
        calls.append(smiles)
        return [{"steps": [], "solved": True}]

    written = run_batch(
        ["CCO", "not_a_real_smiles!!"],
        str(tmp_path),
        timestamp="ts",
        solve=fake_solve,
    )

    assert calls == ["CCO"], "the broken target must never reach the solver"
    bad = written["not_a_real_smiles!!"]
    assert len(bad) == 1 and bad[0].endswith("error.json")
    assert "does not parse" in json.loads(open(bad[0]).read())["error"]
    assert written["CCO"][0].endswith("pathway_1.json")
