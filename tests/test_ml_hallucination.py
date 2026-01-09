"""
Tests for ML-based hallucination detection using XGBoost model.

Tests cover:
- Model loading and caching
- Prediction functionality
- ML hallucination checker integration
- Error handling and fallback behavior
- API endpoint integration
"""

import unittest
import os
import json
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock, Mock
import numpy as np

# Try to import ML hallucination module, skip tests if dependencies missing
try:
    from src.utils.ml_hallucination import (
        load_ml_model,
        predict_hallucination_ml,
        ml_hallucination_checker
    )
    ML_HALLUCINATION_AVAILABLE = True
except ImportError as e:
    ML_HALLUCINATION_AVAILABLE = False
    # Create dummy functions to prevent import errors
    def load_ml_model(*args, **kwargs):
        return False, "Dependencies not available"
    def predict_hallucination_ml(*args, **kwargs):
        return {'is_hallucination': None, 'probability': None, 'method': 'ml_model', 'error': 'Dependencies not available'}
    def ml_hallucination_checker(*args, **kwargs):
        return 500, []
from src.utils.hallucination_checks import hallucination_checker
from src.main import main


@unittest.skipIf(not ML_HALLUCINATION_AVAILABLE, "ML hallucination dependencies not available")
class TestMlHallucinationModelLoading(unittest.TestCase):
    """Test ML model loading functionality."""

    def setUp(self):
        """Set up test fixtures."""
        # Get root directory
        self.root_dir = Path(__file__).parent.parent
        self.models_dir = self.root_dir / "models"
        
        # Store original model files if they exist
        self.original_model_exists = (self.models_dir / "xgboost_hallucination_model.pkl").exists()
        self.original_featurizer_exists = (self.models_dir / "xgboost_featurizer.pkl").exists()
        self.original_metadata_exists = (self.models_dir / "xgboost_metadata.json").exists()

    def tearDown(self):
        """Clean up after tests."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def test_load_model_success(self):
        """Test successful model loading when all files exist."""
        # Only run if model files actually exist
        if not (self.original_model_exists and 
                self.original_featurizer_exists and 
                self.original_metadata_exists):
            self.skipTest("Model files not found. Run training script first.")
        
        success, error = load_ml_model(force_reload=True)
        self.assertTrue(success, f"Model loading failed: {error}")
        self.assertIsNone(error)

    def test_load_model_missing_model_file(self):
        """Test model loading fails gracefully when model file is missing."""
        # Create temporary directory structure
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_models_dir = Path(tmpdir) / "models"
            tmp_models_dir.mkdir()
            
            # Create only featurizer and metadata, not model
            with patch('src.utils.ml_hallucination.Path') as mock_path:
                mock_path.return_value.parent.parent.parent = Path(tmpdir)
                
                success, error = load_ml_model(force_reload=True)
                self.assertFalse(success)
                self.assertIsNotNone(error)
                self.assertIn("Model file not found", error)

    def test_load_model_missing_featurizer_file(self):
        """Test model loading fails gracefully when featurizer file is missing."""
        if not self.original_model_exists:
            self.skipTest("Model file not found. Run training script first.")
        
        # Create temporary directory with only model file
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_models_dir = Path(tmpdir) / "models"
            tmp_models_dir.mkdir()
            
            # Copy only model file
            shutil.copy(
                self.models_dir / "xgboost_hallucination_model.pkl",
                tmp_models_dir / "xgboost_hallucination_model.pkl"
            )
            
            with patch('src.utils.ml_hallucination.Path') as mock_path:
                mock_path.return_value.parent.parent.parent = Path(tmpdir)
                
                success, error = load_ml_model(force_reload=True)
                self.assertFalse(success)
                self.assertIsNotNone(error)
                self.assertIn("Featurizer file not found", error)

    def test_load_model_caching(self):
        """Test that model is cached after first load."""
        if not (self.original_model_exists and 
                self.original_featurizer_exists and 
                self.original_metadata_exists):
            self.skipTest("Model files not found. Run training script first.")
        
        # First load
        success1, error1 = load_ml_model(force_reload=True)
        self.assertTrue(success1)
        
        # Second load should use cache (no force_reload)
        success2, error2 = load_ml_model(force_reload=False)
        self.assertTrue(success2)
        
        # Force reload should reload
        success3, error3 = load_ml_model(force_reload=True)
        self.assertTrue(success3)


@unittest.skipIf(not ML_HALLUCINATION_AVAILABLE, "ML hallucination dependencies not available")
class TestMlHallucinationPrediction(unittest.TestCase):
    """Test ML hallucination prediction functionality."""

    def setUp(self):
        """Set up test fixtures."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def tearDown(self):
        """Clean up after tests."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def test_predict_valid_smiles(self):
        """Test prediction with valid SMILES strings."""
        # Only run if model exists
        root_dir = Path(__file__).parent.parent
        models_dir = root_dir / "models"
        if not (models_dir / "xgboost_hallucination_model.pkl").exists():
            self.skipTest("Model file not found. Run training script first.")
        
        # Test with simple valid transformation
        product = "CCO"  # Ethanol
        reactants = "C=C.O"  # Ethylene + Water
        
        result = predict_hallucination_ml(product, reactants)
        
        self.assertIsNotNone(result)
        self.assertEqual(result['method'], 'ml_model')
        self.assertIn('is_hallucination', result)
        self.assertIn('probability', result)
        self.assertIsInstance(result['is_hallucination'], bool)
        self.assertIsInstance(result['probability'], float)
        self.assertGreaterEqual(result['probability'], 0.0)
        self.assertLessEqual(result['probability'], 1.0)

    def test_predict_invalid_product_smiles(self):
        """Test prediction handles invalid product SMILES."""
        result = predict_hallucination_ml("invalid_smiles", "CCO")
        
        self.assertIsNotNone(result)
        self.assertEqual(result['method'], 'ml_model')
        self.assertTrue(result['is_hallucination'])  # Invalid treated as hallucination
        self.assertEqual(result['probability'], 1.0)
        self.assertIsNotNone(result.get('error'))

    def test_predict_invalid_reactant_smiles(self):
        """Test prediction handles invalid reactant SMILES."""
        result = predict_hallucination_ml("CCO", "invalid_smiles")
        
        self.assertIsNotNone(result)
        self.assertEqual(result['method'], 'ml_model')
        self.assertTrue(result['is_hallucination'])  # Invalid treated as hallucination
        self.assertEqual(result['probability'], 1.0)
        self.assertIsNotNone(result.get('error'))

    def test_predict_identical_molecules(self):
        """Test prediction with identical product and reactants (should be low hallucination)."""
        root_dir = Path(__file__).parent.parent
        models_dir = root_dir / "models"
        if not (models_dir / "xgboost_hallucination_model.pkl").exists():
            self.skipTest("Model file not found. Run training script first.")
        
        smiles = "CCO"  # Ethanol
        result = predict_hallucination_ml(smiles, smiles)
        
        self.assertIsNotNone(result)
        self.assertEqual(result['method'], 'ml_model')
        # Identical molecules should have low hallucination probability
        # (though exact value depends on model training)

    def test_predict_missing_model(self):
        """Test prediction falls back gracefully when model is missing."""
        # Mock load_ml_model to return failure
        with patch('src.utils.ml_hallucination.load_ml_model') as mock_load:
            mock_load.return_value = (False, "Model file not found")
            
            result = predict_hallucination_ml("CCO", "C=C")
            
            self.assertIsNotNone(result)
            self.assertEqual(result['method'], 'ml_model')
            self.assertIsNone(result['is_hallucination'])
            self.assertIsNone(result['probability'])
            self.assertIsNotNone(result.get('error'))


@unittest.skipIf(not ML_HALLUCINATION_AVAILABLE, "ML hallucination dependencies not available")
class TestMlHallucinationChecker(unittest.TestCase):
    """Test ML hallucination checker wrapper function."""

    def setUp(self):
        """Set up test fixtures."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def tearDown(self):
        """Clean up after tests."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def test_ml_checker_valid_pathways(self):
        """Test ML checker with valid pathways."""
        root_dir = Path(__file__).parent.parent
        models_dir = root_dir / "models"
        if not (models_dir / "xgboost_hallucination_model.pkl").exists():
            self.skipTest("Model file not found. Run training script first.")
        
        product = "CCO"  # Ethanol
        res_smiles = [["C=C", "O"]]  # Ethylene + Water
        
        status_code, valid_pathways = ml_hallucination_checker(product, res_smiles)
        
        self.assertEqual(status_code, 200)
        self.assertIsInstance(valid_pathways, list)
        # Valid pathways should be filtered based on ML predictions

    def test_ml_checker_invalid_smiles(self):
        """Test ML checker handles invalid SMILES."""
        product = "CCO"
        res_smiles = [["invalid_smiles"]]
        
        status_code, valid_pathways = ml_hallucination_checker(product, res_smiles)
        
        self.assertEqual(status_code, 200)
        self.assertIsInstance(valid_pathways, list)
        # Invalid SMILES should be filtered out

    def test_ml_checker_empty_list(self):
        """Test ML checker with empty pathway list."""
        product = "CCO"
        res_smiles = []
        
        status_code, valid_pathways = ml_hallucination_checker(product, res_smiles)
        
        self.assertEqual(status_code, 200)
        self.assertEqual(valid_pathways, [])

    def test_ml_checker_multiple_pathways(self):
        """Test ML checker with multiple pathways."""
        root_dir = Path(__file__).parent.parent
        models_dir = root_dir / "models"
        if not (models_dir / "xgboost_hallucination_model.pkl").exists():
            self.skipTest("Model file not found. Run training script first.")
        
        product = "CCO"
        res_smiles = [
            ["C=C", "O"],  # Pathway 1
            ["CC", "O"]    # Pathway 2
        ]
        
        status_code, valid_pathways = ml_hallucination_checker(product, res_smiles)
        
        self.assertEqual(status_code, 200)
        self.assertIsInstance(valid_pathways, list)
        self.assertLessEqual(len(valid_pathways), len(res_smiles))


@unittest.skipIf(not ML_HALLUCINATION_AVAILABLE, "ML hallucination dependencies not available")
class TestHallucinationCheckerIntegration(unittest.TestCase):
    """Test integration of ML model with hallucination_checker."""

    def setUp(self):
        """Set up test fixtures."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def tearDown(self):
        """Clean up after tests."""
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def test_hallucination_checker_rule_based_default(self):
        """Test hallucination_checker defaults to rule-based."""
        product = "CCO"
        res_smiles = [["C=C", "O"]]
        
        status_code, valid_pathways = hallucination_checker(product, res_smiles, use_ml_model=False)
        
        self.assertEqual(status_code, 200)
        self.assertIsInstance(valid_pathways, list)

    def test_hallucination_checker_ml_model(self):
        """Test hallucination_checker with ML model option."""
        root_dir = Path(__file__).parent.parent
        models_dir = root_dir / "models"
        if not (models_dir / "xgboost_hallucination_model.pkl").exists():
            self.skipTest("Model file not found. Run training script first.")
        
        product = "CCO"
        res_smiles = [["C=C", "O"]]
        
        status_code, valid_pathways = hallucination_checker(product, res_smiles, use_ml_model=True)
        
        self.assertEqual(status_code, 200)
        self.assertIsInstance(valid_pathways, list)

    def test_hallucination_checker_ml_fallback(self):
        """Test hallucination_checker falls back to rule-based on ML error."""
        # Mock ml_hallucination_checker to raise exception
        with patch('src.utils.hallucination_checks.ml_hallucination_checker') as mock_ml:
            mock_ml.side_effect = Exception("ML model error")
            
            product = "CCO"
            res_smiles = [["C=C", "O"]]
            
            # Should fall back to rule-based
            status_code, valid_pathways = hallucination_checker(product, res_smiles, use_ml_model=True)
            
            self.assertEqual(status_code, 200)
            self.assertIsInstance(valid_pathways, list)


class TestApiIntegration(unittest.TestCase):
    """Test API endpoint integration with ML hallucination method."""

    def setUp(self):
        """Set up test fixtures."""
        import os
        from src.api import app
        
        # Set test API key
        self.original_api_key = os.environ.get('API_KEY')
        os.environ['API_KEY'] = 'test_api_key_for_ml_tests'
        
        self.app = app
        self.app_context = app.app_context()
        self.app_context.push()
        self.client = app.test_client()
        
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    def tearDown(self):
        """Clean up after tests."""
        import os
        self.app_context.pop()
        
        # Restore original API key
        if self.original_api_key:
            os.environ['API_KEY'] = self.original_api_key
        elif 'API_KEY' in os.environ:
            del os.environ['API_KEY']
        
        # Reset module-level cache
        import src.utils.ml_hallucination as ml_module
        ml_module._ml_model = None
        ml_module._featurizer = None
        ml_module._optimal_threshold = None
        ml_module._model_loaded = False

    @patch('src.api.main')
    def test_api_with_rule_based_method(self, mock_main):
        """Test API accepts rule_based hallucination method."""
        mock_main.return_value = {"result": "test"}
        
        response = self.client.post(
            '/api/retrosynthesis',
            headers={'X-API-KEY': 'test_api_key_for_ml_tests'},
            json={
                'smiles': 'CCO',
                'hallucination_check': 'True',
                'hallucination_method': 'rule_based'
            }
        )
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once()
        # Check that hallucination_method was passed
        call_kwargs = mock_main.call_args[1]
        self.assertEqual(call_kwargs.get('hallucination_method'), 'rule_based')

    @patch('src.api.main')
    def test_api_with_ml_model_method(self, mock_main):
        """Test API accepts ml_model hallucination method."""
        mock_main.return_value = {"result": "test"}
        
        response = self.client.post(
            '/api/retrosynthesis',
            headers={'X-API-KEY': 'test_api_key_for_ml_tests'},
            json={
                'smiles': 'CCO',
                'hallucination_check': 'True',
                'hallucination_method': 'ml_model'
            }
        )
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once()
        # Check that hallucination_method was passed
        call_kwargs = mock_main.call_args[1]
        self.assertEqual(call_kwargs.get('hallucination_method'), 'ml_model')

    @patch('src.api.main')
    def test_api_defaults_to_rule_based(self, mock_main):
        """Test API defaults to rule_based when method not specified."""
        mock_main.return_value = {"result": "test"}
        
        response = self.client.post(
            '/api/retrosynthesis',
            headers={'X-API-KEY': 'test_api_key_for_ml_tests'},
            json={
                'smiles': 'CCO',
                'hallucination_check': 'True'
            }
        )
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once()
        # Check that hallucination_method defaults to rule_based
        call_kwargs = mock_main.call_args[1]
        self.assertEqual(call_kwargs.get('hallucination_method'), 'rule_based')

    @patch('src.api.main')
    def test_api_invalid_method_fallback(self, mock_main):
        """Test API falls back to rule_based for invalid method."""
        mock_main.return_value = {"result": "test"}
        
        response = self.client.post(
            '/api/retrosynthesis',
            headers={'X-API-KEY': 'test_api_key_for_ml_tests'},
            json={
                'smiles': 'CCO',
                'hallucination_check': 'True',
                'hallucination_method': 'invalid_method'
            }
        )
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once()
        # Check that invalid method falls back to rule_based
        call_kwargs = mock_main.call_args[1]
        self.assertEqual(call_kwargs.get('hallucination_method'), 'rule_based')


if __name__ == '__main__':
    unittest.main()

