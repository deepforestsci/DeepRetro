"""Unit tests for deepretro.utils.utils_molecule.

Tests molecule utilities including SMILES validation, substructure matching,
molecular properties, and ring detection.
"""

from __future__ import annotations

import pytest

from deepretro.utils import utils_molecule

# Common SMILES for testing
ETHANOL = "CCO"
BENZENE = "c1ccccc1"
CYCLOHEXANE = "C1CCCCC1"
METHANE = "C"
INVALID_SMILES = "not_a_smiles!!!"


class TestIsValidSmiles:
    """Tests for ``utils_molecule.is_valid_smiles``."""

    @pytest.mark.parametrize(
        ("smiles", "expected"),
        [
            (ETHANOL, True),
            (BENZENE, True),
            (CYCLOHEXANE, True),
            (METHANE, True),
            ("CCOC(=O)C", True),  # ethyl acetate
        ],
    )
    def test_valid_smiles_returns_true(self, smiles: str,
                                       expected: bool) -> None:
        assert utils_molecule.is_valid_smiles(smiles) is expected

    @pytest.mark.parametrize(
        "smiles",
        [
            INVALID_SMILES,
            "C(C(C(C",  # unbalanced
            "xyz",
        ],
    )
    def test_invalid_smiles_returns_false(self, smiles: str) -> None:
        assert utils_molecule.is_valid_smiles(smiles) is False


class TestSubstructureMatching:
    """Tests for ``utils_molecule.substructure_matching``."""

    def test_query_in_target_returns_one(self) -> None:
        # Benzene is a substructure of ethylbenzene
        result = utils_molecule.substructure_matching("CCc1ccccc1", BENZENE)
        assert result == 1

    def test_query_not_in_target_returns_zero(self) -> None:
        # Methane is not a substructure of ethanol in the HasSubstructMatch sense
        # (ethanol has C-C-O, methane C is a single atom - HasSubstructMatch
        # matches if query is subgraph of target; C is subgraph of CCO)
        # Use a case where query is definitely not a substructure
        result = utils_molecule.substructure_matching(METHANE, BENZENE)
        assert result == 0

    def test_benzene_in_ethyl_benzene(self) -> None:
        result = utils_molecule.substructure_matching("CCc1ccccc1", BENZENE)
        assert result == 1

    def test_invalid_target_returns_zero(self) -> None:
        result = utils_molecule.substructure_matching(INVALID_SMILES, BENZENE)
        assert result == 0

    def test_invalid_query_returns_zero(self) -> None:
        result = utils_molecule.substructure_matching(BENZENE, INVALID_SMILES)
        assert result == 0


class TestAreMoleculesSame:
    """Tests for ``utils_molecule.are_molecules_same``."""

    def test_same_smiles_returns_true(self) -> None:
        assert utils_molecule.are_molecules_same(ETHANOL, ETHANOL) is True

    def test_canonically_equivalent_returns_true(self) -> None:
        # Different orderings that canonicalize to same structure
        assert utils_molecule.are_molecules_same("CCO", "OCC") is True

    def test_different_molecules_returns_false(self) -> None:
        assert utils_molecule.are_molecules_same(ETHANOL, BENZENE) is False

    def test_invalid_smiles_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES string provided"):
            utils_molecule.are_molecules_same(INVALID_SMILES, ETHANOL)

    def test_invalid_second_smiles_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES string provided"):
            utils_molecule.are_molecules_same(ETHANOL, INVALID_SMILES)


class TestCalcMolWt:
    """Tests for ``utils_molecule.calc_mol_wt``."""

    def test_ethanol_has_positive_weight(self) -> None:
        result = utils_molecule.calc_mol_wt(ETHANOL)
        assert result > 0
        assert isinstance(result, float)

    def test_methane_weight(self) -> None:
        result = utils_molecule.calc_mol_wt(METHANE)
        assert abs(result - 16.04) < 0.1

    def test_invalid_smiles_returns_zero(self) -> None:
        result = utils_molecule.calc_mol_wt(INVALID_SMILES)
        assert result == 0.0


class TestCalcChemicalFormula:
    """Tests for ``utils_molecule.calc_chemical_formula``."""

    def test_ethanol_formula(self) -> None:
        result = utils_molecule.calc_chemical_formula(ETHANOL)
        assert result == "C2H6O"

    def test_methane_formula(self) -> None:
        result = utils_molecule.calc_chemical_formula(METHANE)
        assert result == "CH4"

    def test_invalid_smiles_returns_na(self) -> None:
        result = utils_molecule.calc_chemical_formula(INVALID_SMILES)
        assert result == "N/A"


class TestComputeFingerprint:
    """Tests for ``utils_molecule.compute_fingerprint``."""

    def test_valid_smiles_returns_list(self) -> None:
        result = utils_molecule.compute_fingerprint(ETHANOL)
        assert isinstance(result, list)
        assert len(result) == 2048  # default nBits
        assert all(b in (0, 1) for b in result)

    def test_custom_radius_and_nbits(self) -> None:
        result = utils_molecule.compute_fingerprint(ETHANOL,
                                                    radius=3,
                                                    nBits=1024)
        assert isinstance(result, list)
        assert len(result) == 1024

    def test_invalid_smiles_returns_none(self) -> None:
        result = utils_molecule.compute_fingerprint(INVALID_SMILES)
        assert result is None


class TestValidityCheck:
    """Tests for ``utils_molecule.validity_check``."""

    def test_valid_pathways_preserved(self) -> None:
        molecule = BENZENE
        res_molecules = [["CC(=O)O"]]
        res_explanations = ["test explanation"]
        res_confidence = [0.8]
        pathways, explanations, confidence = utils_molecule.validity_check(
            molecule, res_molecules, res_explanations, res_confidence)
        assert pathways == [["CC(=O)O"]]
        assert len(explanations) == len(pathways)
        assert len(confidence) == len(pathways)

    def test_invalid_smiles_filtered_out(self) -> None:
        molecule = BENZENE
        res_molecules = [[INVALID_SMILES]]
        res_explanations = ["bad pathway"]
        res_confidence = [0.5]
        pathways, explanations, confidence = utils_molecule.validity_check(
            molecule, res_molecules, res_explanations, res_confidence)
        assert len(pathways) == 0
        assert len(explanations) == 0
        assert len(confidence) == 0

    def test_same_as_target_filtered_out(self) -> None:
        molecule = ETHANOL
        res_molecules = [[ETHANOL]]  # same as target
        res_explanations = ["same molecule"]
        res_confidence = [0.9]
        pathways, explanations, confidence = utils_molecule.validity_check(
            molecule, res_molecules, res_explanations, res_confidence)
        assert len(pathways) == 0
        assert explanations == []
        assert confidence == []

    def test_substructure_of_target_filtered_out(self) -> None:
        molecule = ETHANOL
        res_molecules = [["CC"]]
        res_explanations = ["ethyl fragment"]
        res_confidence = [0.6]
        pathways, explanations, confidence = utils_molecule.validity_check(
            molecule, res_molecules, res_explanations, res_confidence)
        assert pathways == []
        assert explanations == []
        assert confidence == []

    def test_single_smiles_string_instead_of_list(self) -> None:
        molecule = BENZENE
        res_molecules = ["CC(=O)O"
                         ]  # acetic acid - valid, not same, not substructure
        res_explanations = ["test"]
        res_confidence = [0.7]
        pathways, explanations, confidence = utils_molecule.validity_check(
            molecule, res_molecules, res_explanations, res_confidence)
        assert len(pathways) == 1
        assert pathways[0] == ["CC(=O)O"]
        assert explanations[0] == "test"
        assert confidence[0] == 0.7


class TestDetectSevenMemberRings:
    """Tests for ``utils_molecule.detect_seven_member_rings``."""

    def test_cycloheptane_has_seven_member_ring(self) -> None:
        # C1CCCCCC1 is cycloheptane
        result = utils_molecule.detect_seven_member_rings("C1CCCCCC1")
        assert result is True

    def test_cyclohexane_no_seven_member_ring(self) -> None:
        result = utils_molecule.detect_seven_member_rings(CYCLOHEXANE)
        assert result is False

    def test_ethanol_no_rings(self) -> None:
        result = utils_molecule.detect_seven_member_rings(ETHANOL)
        assert result is False

    def test_invalid_smiles_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES string provided"):
            utils_molecule.detect_seven_member_rings(INVALID_SMILES)


class TestDetectEightMemberRings:
    """Tests for ``utils_molecule.detect_eight_member_rings``."""

    def test_cyclooctane_has_eight_member_ring(self) -> None:
        result = utils_molecule.detect_eight_member_rings("C1CCCCCCC1")
        assert result is True

    def test_cycloheptane_no_eight_member_ring(self) -> None:
        result = utils_molecule.detect_eight_member_rings("C1CCCCCC1")
        assert result is False

    def test_cyclohexane_no_eight_member_ring(self) -> None:
        result = utils_molecule.detect_eight_member_rings(CYCLOHEXANE)
        assert result is False

    def test_invalid_smiles_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid SMILES string provided"):
            utils_molecule.detect_eight_member_rings(INVALID_SMILES)
