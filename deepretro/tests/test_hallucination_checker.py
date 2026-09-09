"""Unit tests for the heuristic hallucination checker."""

import pytest

from deepretro.algorithms.hallucination_checker import (
    calculate_hallucination_score,
    hallucination_compare_molecules,
    interpret_score,
)

# Real reaction from the DeepRetro dataset:
# pyrazole bromide + cyclohexanone → α-hydroxycyclohexyl pyrazole
PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"

# Simple valid pair: benzene → methoxybenzene (adding -OC)
BENZENE = "c1ccccc1"
METHOXYBENZENE = "c1ccccc1OC"

# Hallucinated pair: atom counts wildly different
ETHANE = "CC"
LARGE_MOLECULE = "c1ccccc1CCC(=O)C"

INVALID_SMILES = "not_a_smiles!!!"


# hallucination_compare_molecules


def test_compare_valid_molecules_returns_both_valid():
    res = hallucination_compare_molecules(BENZENE, METHOXYBENZENE)
    assert res["valid_reactant"] is True
    assert res["valid_product"] is True


def test_compare_invalid_reactant():
    res = hallucination_compare_molecules(INVALID_SMILES, BENZENE)
    assert res["valid_reactant"] is False
    assert "Invalid reactant SMILES string" in res["detected_issues"]


def test_compare_invalid_product():
    res = hallucination_compare_molecules(BENZENE, INVALID_SMILES)
    assert res["valid_product"] is False
    assert "Invalid product SMILES string" in res["detected_issues"]


def test_compare_identical_molecules_no_issues():
    res = hallucination_compare_molecules(BENZENE, BENZENE)
    assert res["atom_count_consistent"] is True
    assert len(res["detected_issues"]) == 0


def test_compare_detects_atom_count_mismatch():
    res = hallucination_compare_molecules(ETHANE, LARGE_MOLECULE)
    assert res["atom_count_consistent"] is False
    assert any("Atom count mismatch" in i for i in res["detected_issues"])


def test_compare_detects_ring_change():
    # benzene (1 ring) vs ethane (0 rings)
    res = hallucination_compare_molecules(BENZENE, ETHANE)
    assert len(res["ring_size_changes"]) > 0


def test_compare_real_reaction_has_issues():
    # real reaction pair — atoms differ so issues are expected
    res = hallucination_compare_molecules(PYRAZOLE_BROMIDE_KETONE, PYRAZOLE_ADDUCT)
    assert res["valid_reactant"] is True
    assert res["valid_product"] is True


# calculate_hallucination_score


def test_score_returns_dict_with_required_keys():
    result = calculate_hallucination_score(BENZENE, METHOXYBENZENE)
    assert "score" in result
    assert "severity" in result
    assert "message" in result


def test_identical_molecules_are_degenerate_not_clean():
    """Reactants identical to the product decompose nothing.

    This used to score 100/low, so the pipeline would accept a node proposing
    the molecule as its own precursor. It is now rejected as degenerate. No row
    in either label set is affected -- the check is purely defensive.
    """
    result = calculate_hallucination_score("c1ccccc1", "c1ccccc1")
    assert result["score"] == 0
    assert result["severity"] == "critical"
    assert "decompose" in result["message"].lower()


def test_score_invalid_smiles_is_zero():
    result = calculate_hallucination_score(INVALID_SMILES, BENZENE)
    assert result["score"] == 0
    assert result["severity"] == "critical"


def test_score_major_mismatch_is_low():
    result = calculate_hallucination_score(ETHANE, LARGE_MOLECULE)
    assert result["score"] < 80


def test_score_range_0_to_100():
    result = calculate_hallucination_score(BENZENE, METHOXYBENZENE)
    assert 0 <= result["score"] <= 100


def test_score_severity_values():
    for reactant, product in [
        (BENZENE, BENZENE),
        (BENZENE, METHOXYBENZENE),
        (ETHANE, LARGE_MOLECULE),
        (INVALID_SMILES, BENZENE),
    ]:
        result = calculate_hallucination_score(reactant, product)
        assert result["severity"] in ("low", "medium", "high", "critical")


# interpret_score


def test_interpret_score_high():
    msg = interpret_score(95)
    assert msg == "Minimal or no structural inconsistencies detected"


def test_interpret_score_zero():
    msg = interpret_score(0)
    assert "unsupported transformation" in msg.lower()


def test_interpret_score_boundary():
    # The fixed legacy display threshold includes a score of exactly 90.
    msg = interpret_score(90)
    assert msg.startswith("Minimal or no structural")


def test_interpret_score_mid():
    msg = interpret_score(50)
    assert "review recommended" in msg.lower()


# --- Golden regression lock (captured before the compare/score split) ---------

import json  # noqa: E402
from collections import Counter  # noqa: E402
from pathlib import Path  # noqa: E402

_GOLDEN = json.loads(
    (Path(__file__).parent / "fixtures" / "hallucination_golden.json").read_text()
)


@pytest.mark.parametrize("case", _GOLDEN, ids=lambda c: c["product"][:20])
def test_golden_regression_lock(case):
    """Default-weight output matches the reviewed v2 fixtures."""
    result = calculate_hallucination_score(case["reactants"], case["product"])
    # The exact key set matters: the invalid-SMILES branch returns only 3 keys.
    assert sorted(result.keys()) == case["keys"]
    assert result["score"] == case["score"]
    assert result["severity"] == case["severity"]
    assert result["message"] == case["message"]
    # penalties order is not deterministic pre-refactor: compare as a multiset.
    assert Counter(result.get("penalties", [])) == Counter(case["penalties"])


def test_output_is_hash_seed_independent():
    """After the sorted() fix, even list ordering must be stable."""
    import os
    import subprocess
    import sys

    snippet = (
        "from deepretro.algorithms.hallucination_checker import "
        "calculate_hallucination_score as f;"
        "r=f('CCOCCN','ClC(=O)c1ccccc1Br');"
        "print(r['score'], r['severity'], '|'.join(r['penalties']))"
    )
    outs = set()
    for seed in ("0", "1", "2", "3", "4"):
        # Inherit the environment and override only the seed: a bare env would
        # strip the venv and fail to import deepretro.
        env = dict(os.environ, PYTHONHASHSEED=seed)
        proc = subprocess.run(
            [sys.executable, "-c", snippet],
            capture_output=True,
            text=True,
            check=True,
            env=env,
        )
        outs.add(proc.stdout.strip())
    assert len(outs) == 1, f"output varies with PYTHONHASHSEED: {outs}"


# --- Carbon radicals are truncated SMILES, not chemistry --------------------


def test_carbon_radical_is_treated_as_unassessable():
    """Measured: fires on 8 clean human rows, all 8 hallucinated.

    Five were accepted at scores 40/85/90/95/95 before this, because the atom
    and ring counts of a truncated SMILES look fine.
    """
    result = calculate_hallucination_score("c1ccccc1[CH2]", "c1ccccc1C")
    assert result["score"] == 0
    assert result["severity"] == "critical"
    assert "radical" in result["message"].lower()


@pytest.mark.parametrize(
    "smiles",
    [
        "CC1(C)CCCC(C)(C)N1[O]",  # TEMPO
        "CC(C)(C)[O]",  # t-butoxyl
        "[N+](=O)[O-]",  # nitro
        "O=O",  # dioxygen
    ],
)
def test_heteroatom_radicals_are_not_flagged(smiles):
    """Restricting to carbon is what keeps this at zero false positives."""
    comparison = hallucination_compare_molecules(smiles, smiles)
    assert comparison.get("carbon_radical") is False


def test_carbon_radical_flag_is_symmetric_over_reactant_and_product():
    for reactant, product in (
        ("c1ccccc1[CH2]", "c1ccccc1C"),
        ("c1ccccc1C", "c1ccccc1[CH2]"),
    ):
        assert (
            hallucination_compare_molecules(reactant, product)["carbon_radical"] is True
        )


@pytest.mark.parametrize(
    "name,reactant,product",
    [
        ("toluene nitration", "Cc1ccccc1", "Cc1ccc([N+](=O)[O-])cc1"),
        ("phenol bromination", "Oc1ccccc1", "Oc1ccc(Br)cc1"),
        ("Friedel-Crafts acylation", "c1ccccc1.CC(=O)Cl", "CC(=O)c1ccccc1"),
    ],
)
def test_electrophilic_aromatic_substitution_is_not_rejected(name, reactant, product):
    """Adding a substituent relabels the others, which is not a rearrangement.

    Positions are relative (ortho/meta/para), so nitrating toluene "moves" the
    methyl and the check fired on every C-H functionalisation of an arene.
    Nitration scored 15 and bromination 30, both below the accept threshold.
    """
    assert calculate_hallucination_score(reactant, product)["score"] >= 40, name


def test_a_genuine_position_change_is_still_caught():
    """The check must not be defanged: same substituent count, moved."""
    result = calculate_hallucination_score("Cc1ccc(O)cc1", "Cc1cccc(O)c1")
    assert result["score"] < 40
    assert result["severity"] in ("high", "critical")


def test_score_does_not_depend_on_reactant_fragment_order():
    """Ring pairing took matching_rings[0], so the verdict was order-dependent."""
    a = calculate_hallucination_score("Cc1ccc(O)cc1.Cc1ccccc1Cl", "Cc1ccc(O)cc1")[
        "score"
    ]
    b = calculate_hallucination_score("Cc1ccccc1Cl.Cc1ccc(O)cc1", "Cc1ccc(O)cc1")[
        "score"
    ]
    assert a == b
