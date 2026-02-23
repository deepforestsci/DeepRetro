"""Unit tests for the heuristic hallucination checker."""

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


def test_score_identical_molecules_is_high():
    result = calculate_hallucination_score(BENZENE, BENZENE)
    assert result["score"] == 100
    assert result["severity"] == "low"


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
    assert "reliable" in msg.lower()


def test_interpret_score_zero():
    msg = interpret_score(0)
    assert "hallucination" in msg.lower()


def test_interpret_score_boundary():
    # score of exactly 90 should be "Highly reliable"
    msg = interpret_score(90)
    assert msg.startswith("Highly reliable")


def test_interpret_score_mid():
    msg = interpret_score(50)
    assert "questionable" in msg.lower()

