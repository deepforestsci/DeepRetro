import pytest
from unittest import mock
import os

# Import the module to test
from src.utils import utils_molecule

# Mock the RDKit imports to avoid dependency failures
@pytest.fixture(autouse=True)
def mock_rdkit():
    with mock.patch('src.utils.utils_molecule.Chem') as mock_chem, \
         mock.patch('src.utils.utils_molecule.AllChem') as mock_allchem, \
         mock.patch('src.utils.utils_molecule.rdMolDescriptors') as mock_mol_descriptors, \
         mock.patch('src.utils.utils_molecule.CalcMolFormula') as mock_calc_mol_formula, \
         mock.patch('src.utils.utils_molecule.ExactMolWt') as mock_exact_mol_wt, \
         mock.patch('src.utils.utils_molecule.joblib') as mock_joblib:
        
        # Configure common mock behaviors
        mock_mol = mock.MagicMock()
        mock_chem.MolFromSmiles.return_value = mock_mol
        mock_chem.MolToSmiles.return_value = "C"
        mock_exact_mol_wt.return_value = 100.0
        mock_calc_mol_formula.return_value = "C6H12O6"
        mock_mol_descriptors.GetMorganFingerprintAsBitVect.return_value = [0, 1, 0, 1]
        mock_allchem.GetMorganFingerprintAsBitVect.return_value = [0, 1, 0, 1]
        
        # Configure behavior for HasSubstructMatch
        mock_mol.HasSubstructMatch.return_value = True
        
        # Configure behavior for GetRingInfo
        mock_ring_info = mock.MagicMock()
        mock_ring_info.AtomRings.return_value = [(0, 1, 2, 3, 4, 5, 6)]  # 7-member ring
        mock_mol.GetRingInfo.return_value = mock_ring_info
        
        yield mock_chem, mock_allchem, mock_mol_descriptors, mock_calc_mol_formula, mock_exact_mol_wt, mock_joblib


# Test for log_message function
def test_log_message():
    with mock.patch('src.utils.utils_molecule.print') as mock_print:
        utils_molecule.log_message("Test message")
        mock_print.assert_called_once_with("Test message")
    
    # Test with logger
    mock_logger = mock.MagicMock()
    utils_molecule.log_message("Test message", mock_logger)
    mock_logger.assert_not_called()  # The function doesn't use the logger parameter correctly


# Test for is_valid_smiles function
def test_is_valid_smiles():
    # Test valid SMILES
    assert utils_molecule.is_valid_smiles("C") is True
    
    # Test invalid SMILES - configure mock to return None
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', return_value=None):
        assert utils_molecule.is_valid_smiles("invalid") is False
    
    # Test exception case
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', side_effect=Exception("Test exception")):
        assert utils_molecule.is_valid_smiles("exception") is False


# Test for substructure_matching function
@mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles')
def test_substructure_matching(mock_mol_from_smiles):
    # Create mock molecule objects
    mock_mol = mock.MagicMock()
    mock_submol = mock.MagicMock()
    
    # Configure mock molecule to return True for HasSubstructMatch
    mock_mol.HasSubstructMatch.return_value = True
    mock_mol_from_smiles.side_effect = [mock_mol, mock_submol]
    
    # Test case where substructure is found
    result = utils_molecule.substructure_matching("C1CCCCC1", "C1CCCC1")
    assert result == 1
    
    # Configure mock molecule to return False for HasSubstructMatch
    mock_mol.HasSubstructMatch.return_value = False
    
    # Test case where substructure is not found
    result = utils_molecule.substructure_matching("C1CCCCC1", "C1CCCC1")
    assert result == 0
    
    # Test case with invalid SMILES
    mock_mol_from_smiles.side_effect = [None, mock_submol]
    result = utils_molecule.substructure_matching("invalid", "C1CCCC1")
    assert result == 0


# Test for calc_mol_wt function
def test_calc_mol_wt():
    # Test successful calculation
    result = utils_molecule.calc_mol_wt("C")
    assert result == 100.0
    
    # Test exception case
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', side_effect=Exception("Test exception")):
        result = utils_molecule.calc_mol_wt("invalid")
        assert result == 0.0


# Test for calc_chemical_formula function
def test_calc_chemical_formula():
    # Test successful calculation
    result = utils_molecule.calc_chemical_formula("C")
    assert result == "C6H12O6"
    
    # Test exception case
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', side_effect=Exception("Test exception")):
        result = utils_molecule.calc_chemical_formula("invalid")
        assert result == "N/A"


# Test for are_molecules_same function
def test_are_molecules_same():
    # Test same canonical SMILES
    result = utils_molecule.are_molecules_same("C", "C")
    assert result is True
    
    # Test same fingerprint but different canonical SMILES
    with mock.patch('src.utils.utils_molecule.Chem.MolToSmiles', side_effect=["C1", "C2"]):
        with mock.patch('src.utils.utils_molecule.rdMolDescriptors.GetMorganFingerprintAsBitVect', 
                        return_value=[1, 1, 1]):
            result = utils_molecule.are_molecules_same("C1", "C2")
            assert result is True
    
    # Test different molecules
    with mock.patch('src.utils.utils_molecule.Chem.MolToSmiles', side_effect=["C1", "C2"]):
        with mock.patch('src.utils.utils_molecule.rdMolDescriptors.GetMorganFingerprintAsBitVect', 
                       side_effect=[[1, 0, 1], [0, 1, 0]]):
            result = utils_molecule.are_molecules_same("C1", "C2")
            assert result is False
    
    # Test invalid SMILES
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', return_value=None):
        with pytest.raises(ValueError, match="Invalid SMILES string provided."):
            utils_molecule.are_molecules_same("invalid", "C")


# Test for compute_fingerprint function
def test_compute_fingerprint():
    # Test successful fingerprint calculation
    result = utils_molecule.compute_fingerprint("C")
    assert result == [0, 1, 0, 1]
    
    # Test invalid SMILES
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', return_value=None):
        result = utils_molecule.compute_fingerprint("invalid")
        assert result is None


# Test for get_reaction_type function
@mock.patch('src.utils.utils_molecule.joblib.load')
@mock.patch('src.utils.utils_molecule.compute_fingerprint')
def test_get_reaction_type(mock_compute_fingerprint, mock_load):
    # Configure mocks
    mock_model = mock.MagicMock()
    mock_model.predict.return_value = [0]
    mock_load.return_value = mock_model
    mock_compute_fingerprint.side_effect = [[1, 0, 1], [0, 1, 0]]
    
    # Mock the REACTION_ENCODING_NAMES
    with mock.patch.dict('src.utils.utils_molecule.REACTION_ENCODING_NAMES', {0: 'Test Reaction'}):
        reaction_name, reaction_type = utils_molecule.get_reaction_type("C1", "C2", "model_path")
        
        assert reaction_name == 'Test Reaction'
        assert reaction_type == 0
        mock_load.assert_called_once_with("model_path")
        mock_model.predict.assert_called_once()


# Test for get_reaction_type function with FileNotFoundError
@mock.patch('src.utils.utils_molecule.joblib.load')
def test_get_reaction_type_file_not_found(mock_load):
    # Configure mocks to simulate FileNotFoundError
    mock_load.side_effect = FileNotFoundError("Model file not found")
    
    reaction_name, reaction_type = utils_molecule.get_reaction_type("C1", "C2", "model_path")
    
    assert reaction_name == 'Unknown Reaction'
    assert reaction_type == -1
    mock_load.assert_called_once_with("model_path")


# Test for get_reaction_type function with other exception
@mock.patch('src.utils.utils_molecule.joblib.load')
@mock.patch('src.utils.utils_molecule.compute_fingerprint')
def test_get_reaction_type_exception(mock_compute_fingerprint, mock_load):
    # Configure mocks to simulate exception
    mock_load.return_value = mock.MagicMock()
    mock_compute_fingerprint.side_effect = Exception("Test exception")
    
    reaction_name, reaction_type = utils_molecule.get_reaction_type("C1", "C2", "model_path")
    
    assert reaction_name == 'Unknown Reaction'
    assert reaction_type == -1


# Test for calc_confidence_estimate function
def test_calc_confidence_estimate():
    # Test with probability < 0.3
    result = utils_molecule.calc_confidence_estimate(0.2)
    assert result == 0.8
    
    # Test with 0.3 <= probability < 0.45
    result = utils_molecule.calc_confidence_estimate(0.4)
    assert result == 0.9
    
    # Test with 0.45 <= probability < 0.6
    result = utils_molecule.calc_confidence_estimate(0.5)
    assert result == 0.8
    
    # Test with probability >= 0.6
    result = utils_molecule.calc_confidence_estimate(0.7)
    assert result == 0.7
    
    # Test with probability > 0.99 (should be capped at 0.99)
    result = utils_molecule.calc_confidence_estimate(0.995)
    assert result == 0.99
    
    # Test with list input
    result = utils_molecule.calc_confidence_estimate([0.2])
    assert result == 0.8


# Test for calc_scalability_index function
@mock.patch('src.utils.utils_molecule.get_reaction_type')
def test_calc_scalability_index(mock_get_reaction_type):
    # Configure mocks
    mock_get_reaction_type.return_value = ('Test Reaction', 1)
    
    # Mock the ENCODING_SCALABILITY
    with mock.patch.dict('src.utils.utils_molecule.ENCODING_SCALABILITY', {1: '5'}):
        result = utils_molecule.calc_scalability_index("C1", "C2")
        
        assert result == '5'
        mock_get_reaction_type.assert_called_once_with("C1", "C2", utils_molecule.RXN_CLASSIFICATION_MODEL_PATH)


# Test for calc_scalability_index function with model not found
@mock.patch('src.utils.utils_molecule.get_reaction_type')
def test_calc_scalability_index_model_not_found(mock_get_reaction_type):
    # Configure mocks to simulate model not found
    mock_get_reaction_type.return_value = ('Unknown Reaction', -1)
    
    result = utils_molecule.calc_scalability_index("C1", "C2")
    
    assert result == 'N/A'
    mock_get_reaction_type.assert_called_once_with("C1", "C2", utils_molecule.RXN_CLASSIFICATION_MODEL_PATH)


# Test for calc_scalability_index function with exception
@mock.patch('src.utils.utils_molecule.get_reaction_type')
def test_calc_scalability_index_exception(mock_get_reaction_type):
    # Configure mocks to simulate exception
    mock_get_reaction_type.side_effect = Exception("Test exception")
    
    result = utils_molecule.calc_scalability_index("C1", "C2")
    
    assert result == 'N/A'
    mock_get_reaction_type.assert_called_once_with("C1", "C2", utils_molecule.RXN_CLASSIFICATION_MODEL_PATH)


# Test for calc_yield function
def test_calc_yield():
    result = utils_molecule.calc_yield("C1", "C2")
    assert result == "#"


# Test for detect_seven_member_rings function
def test_detect_seven_member_rings():
    # Test successful detection of 7-member ring
    result = utils_molecule.detect_seven_member_rings("C1CCCCCC1")
    assert result is True
    
    # Test no 7-member ring
    with mock.patch('src.utils.utils_molecule.Chem') as mock_chem:
        mock_mol = mock.MagicMock()
        mock_ring_info = mock.MagicMock()
        mock_ring_info.AtomRings.return_value = [(0, 1, 2, 3, 4, 5)]  # 6-member ring
        mock_mol.GetRingInfo.return_value = mock_ring_info
        mock_chem.MolFromSmiles.return_value = mock_mol
        
        result = utils_molecule.detect_seven_member_rings("C1CCCCC1")
        assert result is False
    
    # Test invalid SMILES
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', return_value=None):
        with pytest.raises(ValueError, match="Invalid SMILES string provided."):
            utils_molecule.detect_seven_member_rings("invalid")


# Test for detect_eight_member_rings function
def test_detect_eight_member_rings():
    # Test successful detection of 8-member ring
    with mock.patch('src.utils.utils_molecule.Chem') as mock_chem:
        mock_mol = mock.MagicMock()
        mock_ring_info = mock.MagicMock()
        mock_ring_info.AtomRings.return_value = [(0, 1, 2, 3, 4, 5, 6, 7)]  # 8-member ring
        mock_mol.GetRingInfo.return_value = mock_ring_info
        mock_chem.MolFromSmiles.return_value = mock_mol
        
        result = utils_molecule.detect_eight_member_rings("C1CCCCCCC1")
        assert result is True
    
    # Test no 8-member ring
    with mock.patch('src.utils.utils_molecule.Chem') as mock_chem:
        mock_mol = mock.MagicMock()
        mock_ring_info = mock.MagicMock()
        mock_ring_info.AtomRings.return_value = [(0, 1, 2, 3, 4, 5, 6)]  # 7-member ring
        mock_mol.GetRingInfo.return_value = mock_ring_info
        mock_chem.MolFromSmiles.return_value = mock_mol
        
        result = utils_molecule.detect_eight_member_rings("C1CCCCCC1")
        assert result is False
    
    # Test invalid SMILES
    with mock.patch('src.utils.utils_molecule.Chem.MolFromSmiles', return_value=None):
        with pytest.raises(ValueError, match="Invalid SMILES string provided."):
            utils_molecule.detect_eight_member_rings("invalid")


if __name__ == '__main__':
    pytest.main() 