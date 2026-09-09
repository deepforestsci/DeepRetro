"""Regression tests for isotope handling, ranking and structural comparisons."""

import json
from typing import Any

import pytest

from deepretro.algorithms.hallucination_checker import (
    calculate_hallucination_score,
    hallucination_compare_molecules,
)
from deepretro.algorithms.pipeline_checks import hallucination_checker

# --- F1: isotope labelling is not a degenerate step -------------------------


@pytest.mark.parametrize(
    "reactant,product",
    [
        ("[2H]C([2H])=Cc1ccccc1", "C=Cc1ccccc1"),  # deuterated tracer
        ("[13CH3]c1ccccc1", "Cc1ccccc1"),  # 13C internal standard
    ],
)
def test_isotope_labelling_is_not_degenerate(reactant: Any, product: Any) -> None:
    """F1. Normalising isotopes made every labelling step look self-identical.

    Deuterated tracers and 13C standards are real synthetic chemistry; the whole
    class was scoring 0/critical. Degeneracy is now judged before normalisation.
    """
    result = calculate_hallucination_score(reactant, product)
    assert result["score"] >= 40
    assert "Degenerate" not in result["message"]


def test_a_real_degenerate_step_is_still_caught() -> None:
    """F1 must not be fixed by disabling the degeneracy check."""
    result = calculate_hallucination_score("CCO", "CCO")
    assert result["score"] == 0
    assert result["severity"] == "critical"
    assert result["unassessable"] is True


# --- F2: a discarded fragment must not change the verdict -------------------


@pytest.mark.parametrize(
    "reactant,product",
    [
        ("[2H]C([2H])=Cc1ccccc1", "C=Cc1ccccc1"),
        ("[13CH3]c1ccccc1", "Cc1ccccc1"),
    ],
)
def test_an_ignored_spectator_cannot_change_the_score(
    reactant: Any, product: Any
) -> None:
    """F2. Stripping re-parsed the reactants and threw away the normalisation.

    d2-styrene scored 0 alone but 90 with an ethyl acetate appended -- a fragment
    the checker itself then discards. A 90-point swing from irrelevant input.
    """
    alone = calculate_hallucination_score(reactant, product)["score"]
    with_spectator = calculate_hallucination_score(f"{reactant}.CCOC(C)=O", product)[
        "score"
    ]
    assert alone == with_spectator


# --- F5: hydrogen is not a reactant -----------------------------------------


def test_molecular_hydrogen_is_a_spectator() -> None:
    """F5. [H][H] parses but has no heavy atoms, so it reached the counter.

    Hydrogen is invisible to the checker everywhere else, which made catalytic
    hydrogenation the one reaction taxed for the H2 it consumes.
    """
    comparison = hallucination_compare_molecules("C=Cc1ccccc1.[H][H]", "CCc1ccccc1")
    assert comparison["atom_count_deltas"] == {}
    assert (
        calculate_hallucination_score("C=Cc1ccccc1.[H][H]", "CCc1ccccc1")["score"]
        == 100
    )


# --- F6: the demoted tier must be ordered -----------------------------------


def test_penalty_total_orders_candidates_whose_score_has_floored() -> None:
    """F6. `score` floors at 0, where 51% of demoted candidates tied.

    A demoted tier is only useful if it is ordered, and the bottom of the scale
    is exactly where ordering matters most.
    """
    mild = calculate_hallucination_score("CCO", "c1ccc2ccccc2c1")
    severe = calculate_hallucination_score("C", "c1ccc2ccccc2c1CCCCCCCCCC")
    assert mild["score"] == severe["score"] == 0
    assert severe["penalty_total"] > mild["penalty_total"]


def test_penalty_total_is_zero_for_a_clean_step() -> None:
    assert calculate_hallucination_score("CCO", "CC=O")["penalty_total"] == 0


# --- F7: unassessable candidates are dropped, not ranked --------------------


def test_unassessable_candidates_are_never_returned() -> None:
    """F7. Three hard failures all score 0 and were offered as last resorts.

    The target was being handed back as its own precursor, and a truncated
    SMILES alongside it. The checker's own verdict on these is "cannot be
    assessed"; ranking must not override that into "try this last".
    """
    status, ranked = hallucination_checker(
        "CCc1ccccc1",
        [["CCc1ccccc1"], ["c1ccccc1[CH2]"], ["C=Cc1ccccc1"]],
    )
    assert status == 200
    assert ranked == [["C=Cc1ccccc1"]]


def test_a_merely_suspicious_candidate_is_still_demoted_not_dropped() -> None:
    """The third tier must not swallow the demoted tier."""
    status, ranked = hallucination_checker("c1ccc2ccccc2c1", [["CCO"]])
    assert status == 200
    assert ranked == [["CCO"]]


# --- F9: the message must agree with the severity ---------------------------


def test_message_cannot_contradict_severity() -> None:
    """F9. interpret_score hardcoded thresholds while severity used the tunable ones.

    That produced reports reading "score 30, severity low, Likely hallucination".
    """
    from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS

    weights = DEFAULT_WEIGHTS.replace(cut_low=25, cut_medium=22, cut_high=20)
    result = calculate_hallucination_score("CCCCCCCCCCCCCCC", "C", weights=weights)
    assert result["severity"] == "low"
    assert "hallucination" not in result["message"].lower()


# --- F10: one dict shape on every path --------------------------------------


def test_comparison_dict_has_a_stable_shape() -> None:
    """F10. These are written to feature files and JSONL."""
    clean = hallucination_compare_molecules("CCO", "CC=O")
    stripped = hallucination_compare_molecules(
        "O=Cc1ccccc1.C=P(c1ccccc1)(c1ccccc1)c1ccccc1", "C=Cc1ccccc1"
    )
    invalid = hallucination_compare_molecules("not_a_smiles", "CCO")
    assert set(clean) == set(stripped)
    assert set(invalid) <= set(clean)
    for key in ("degenerate", "carbon_radical", "dropped_reagents"):
        assert key in clean and key in invalid


# --- F12: the verdict must survive into the delivered route -----------------


def test_verdict_survives_format_output() -> None:
    """F12. format_output rebuilds the schema and dropped reaction metadata.

    The verdict was recorded on the tree and then discarded before anything
    could read it, so a route leaning on a flagged step was still silent -- the
    exact failure the recording was meant to prevent.
    """
    from deepretro.algorithms.autosolve import AutoSolver, reaction_tree, unsolved_leaf

    tree = reaction_tree(
        "CCc1ccccc1",
        [unsolved_leaf("C=Cc1ccccc1")],
        [0.8],
        {"score": 35, "severity": "high", "flagged": True},
    )
    output = AutoSolver(hallucination_mode="none").parse(tree, solved=True)

    assert "hallucination" in json.dumps(output)
    summary = output["hallucination_summary"]
    assert summary["n_flagged"] == 1
    assert summary["min_score"] == 35
    assert summary["flagged_steps"][0]["product"] == "CCc1ccccc1"


def test_gate_off_routes_report_no_verdicts() -> None:
    from deepretro.algorithms.autosolve import AutoSolver, reaction_tree, unsolved_leaf

    tree = reaction_tree("CC=O", [unsolved_leaf("CCO")], [0.9])
    summary = AutoSolver(hallucination_mode="none").parse(tree, solved=True)[
        "hallucination_summary"
    ]
    assert summary == {
        "n_steps": 0,
        "n_flagged": 0,
        "min_score": None,
        "flagged_steps": [],
    }


# --- F8: the printed penalties must sum to the applied score ----------------


def test_printed_penalties_sum_to_the_applied_penalty_when_capped() -> None:
    """F8. Per-element lines summed to more than the score once cap_atom bound.

    A case printing -25, -30 and -30 had only 30 applied: the explanation
    overstated the penalty by 55 points.
    """
    import re

    from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS

    result = calculate_hallucination_score(
        "CCCCCCNNNNNNOOOOOO", "C", weights=DEFAULT_WEIGHTS.replace(cap_atom=30)
    )
    printed = sum(
        int(m)
        for line in result["penalties"]
        for m in re.findall(r"-(\d+) points", line)
    )
    assert printed == 100 - result["score"]


def test_uncapped_case_still_lists_each_element() -> None:
    result = calculate_hallucination_score("CCO", "CCN")
    assert len(result["penalties"]) == 2
    assert all("capped" not in line for line in result["penalties"])


# --- F13: the verdict must name the gate that produced it -------------------


def test_verdict_names_its_source() -> None:
    """F13. _verdict_for ran the heuristic regardless of the configured mode.

    An ML-gated route was annotated with the heuristic's score, so `flagged`
    meant "the heuristic would have rejected this", not "the gate that actually
    ran did".
    """
    from deepretro.algorithms.autosolve import AutoSolver

    verdict = AutoSolver(hallucination_mode="heuristic")._verdict_for(["CCO"], "CC=O")
    assert verdict["source"] == "heuristic"


# --- F3: stripping must not hide a hallucination ----------------------------


@pytest.mark.parametrize(
    "name,reactants",
    [
        ("benzene + dodecane", "c1ccccc1.CCCCCCCCCCCC"),
        ("benzene + steroid", "c1ccccc1.CC(C)CCCC(C)C1CCC2C1CCC1C2CCC2CCCCC12"),
    ],
)
def test_stripping_cannot_rescue_reactants_too_small_for_the_product(
    name: Any, reactants: Any
) -> None:
    """F3. An invented fragment has the same signature as a spectator reagent.

    Coverage alone deleted exactly the fragments that produce the checker's
    strongest signal: benzene plus a steroid skeleton proposed for ethylbenzene
    scored 0 without stripping and 80 with it. A fragment may now only be
    dropped if what remains can still account for the product.
    """
    assert calculate_hallucination_score(reactants, "CCc1ccccc1")["score"] < 60, name


@pytest.mark.parametrize(
    "name,reactants,product",
    [
        ("Wittig", "O=Cc1ccccc1.C=P(c1ccccc1)(c1ccccc1)c1ccccc1", "C=Cc1ccccc1"),
        (
            "DCC coupling",
            "CC(=O)O.NCc1ccccc1.C(=NC1CCCCC1)=NC1CCCCC1",
            "CC(=O)NCc1ccccc1",
        ),
        ("aldol", "CC=O.CC=O", "CC(O)CC=O"),
        ("Suzuki", "OB(O)c1ccccc1.Brc1ccccc1", "c1ccc(-c2ccccc2)cc1"),
    ],
)
def test_real_reagents_are_still_ignored(
    name: Any, reactants: Any, product: Any
) -> None:
    """The guard must not undo the fix it protects."""
    assert calculate_hallucination_score(reactants, product)["score"] >= 40, name


# --- F11: the position check must work on non-six-membered rings ------------


def test_position_change_is_detected_on_a_five_membered_ring() -> None:
    """F11. determine_ring_position returned a constant "1" for any non-6-ring.

    Every substituent on a pyrrole, furan, imidazole, thiophene or pyrazole got
    the same label, so the comparison could never fire -- w_pos was inert across
    a large slice of medicinal chemistry, and part of its tuning range was
    meaningless.
    """
    # Adjacent -> separated: the two substituents' ring separation changes from
    # 1 to 2, which the check can now see. Before, both were labelled "1".
    moved = calculate_hallucination_score("Cc1ccoc1C", "Cc1cc(C)co1")
    assert moved["score"] < 100
    assert moved["n_position_changes"] if "n_position_changes" in moved else True


def test_five_ring_position_check_only_sees_separation_not_identity() -> None:
    """Honest limit of the generalisation.

    The label is the minimum ring-bond separation to another substituent, which
    generalises ortho/meta/para. On a five-ring, 2,5- and 2,4-disubstitution
    both have minimum separation 2, so this pair is indistinguishable. Detecting
    it would need a canonical numbering, not a distance.
    """
    same_separation = calculate_hallucination_score("Cc1ccc(C)o1", "Cc1cc(C)co1")
    assert same_separation["score"] == 100


def test_six_membered_ring_behaviour_is_unchanged() -> None:
    assert calculate_hallucination_score("Cc1ccc(O)cc1", "Cc1cccc(O)c1")["score"] < 40
    assert (
        calculate_hallucination_score("Cc1ccccc1", "Cc1ccc([N+](=O)[O-])cc1")["score"]
        >= 40
    )
