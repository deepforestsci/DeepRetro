"""
ML-based hallucination detection using trained XGBoost model.

This module provides functions to load and use the trained XGBoost model
for hallucination prediction, as an alternative to the rule-based system.
"""

import pickle
import json
from pathlib import Path
import numpy as np
from typing import Optional, Dict, Tuple

from src.utils.utils_molecule import is_valid_smiles
from src.utils.job_context import logger as context_logger
import os

ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING",
                                    "true").lower() == "false" else True

# Lazy imports - only import when needed
_CircularFingerprint = None
_XGBClassifier = None


def _import_dependencies():
    """Lazy import of deepchem and xgboost dependencies."""
    global _CircularFingerprint, _XGBClassifier
    if _CircularFingerprint is None:
        try:
            from deepchem.feat import CircularFingerprint
            _CircularFingerprint = CircularFingerprint
        except ImportError:
            _CircularFingerprint = False  # Mark as unavailable
    
    if _XGBClassifier is None:
        try:
            from xgboost import XGBClassifier
            _XGBClassifier = XGBClassifier
        except ImportError:
            _XGBClassifier = False  # Mark as unavailable
    
    return _CircularFingerprint is not False and _XGBClassifier is not False


def log_message(message: str, logger=None):
    """Log a message using the context logger if available, otherwise print."""
    if ENABLE_LOGGING and logger:
        logger.info(message)
    else:
        print(message)


# Global variables to cache loaded model (lazy loading)
_ml_model: Optional[object] = None  # XGBClassifier when loaded
_featurizer: Optional[object] = None  # CircularFingerprint when loaded
_optimal_threshold: Optional[float] = None
_model_loaded: bool = False


def load_ml_model(force_reload: bool = False) -> Tuple[bool, Optional[str]]:
    """
    Load the trained XGBoost model, featurizer, and metadata.
    
    Uses lazy loading with caching - model is loaded once and reused.
    
    Parameters
    ----------
    force_reload : bool, optional
        If True, force reload even if model is already loaded (default: False)
    
    Returns
    -------
    Tuple[bool, Optional[str]]
        (success, error_message)
        - success: True if model loaded successfully, False otherwise
        - error_message: Error description if loading failed, None if successful
    """
    global _ml_model, _featurizer, _optimal_threshold, _model_loaded
    
    # Return cached model if already loaded and not forcing reload
    if _model_loaded and not force_reload:
        return True, None
    
    logger = context_logger.get() if ENABLE_LOGGING else None
    
    try:
        # Get root directory (assuming we're in src/utils/)
        root_dir = Path(__file__).parent.parent.parent
        models_dir = root_dir / "models"
        
        model_path = models_dir / "xgboost_hallucination_model.pkl"
        featurizer_path = models_dir / "xgboost_featurizer.pkl"
        metadata_path = models_dir / "xgboost_metadata.json"
        
        # Check if all required files exist
        if not model_path.exists():
            error_msg = f"Model file not found: {model_path}"
            log_message(f"ML Model: {error_msg}", logger)
            return False, error_msg
        
        if not featurizer_path.exists():
            error_msg = f"Featurizer file not found: {featurizer_path}"
            log_message(f"ML Model: {error_msg}", logger)
            return False, error_msg
        
        if not metadata_path.exists():
            error_msg = f"Metadata file not found: {metadata_path}"
            log_message(f"ML Model: {error_msg}", logger)
            return False, error_msg
        
        # Load model
        with open(model_path, 'rb') as f:
            _ml_model = pickle.load(f)
        log_message(f"ML Model: Loaded model from {model_path}", logger)
        
        # Load featurizer (it should be a CircularFingerprint instance)
        with open(featurizer_path, 'rb') as f:
            _featurizer = pickle.load(f)
        # Verify it's the right type
        if not isinstance(_featurizer, _CircularFingerprint):
            return False, f"Featurizer is not a CircularFingerprint instance"
        log_message(f"ML Model: Loaded featurizer from {featurizer_path}", logger)
        
        # Load metadata
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        _optimal_threshold = metadata.get('optimal_threshold', 0.5)
        log_message(f"ML Model: Loaded metadata, optimal threshold = {_optimal_threshold:.4f}", logger)
        
        _model_loaded = True
        return True, None
        
    except Exception as e:
        error_msg = f"Error loading ML model: {str(e)}"
        log_message(f"ML Model: {error_msg}", logger)
        _model_loaded = False
        return False, error_msg


def predict_hallucination_ml(product_smiles: str, reactants_smiles: str) -> Dict:
    """
    Predict hallucination using the ML model.
    
    Parameters
    ----------
    product_smiles : str
        SMILES string of the product molecule
    reactants_smiles : str
        SMILES string of reactants (can be multiple, separated by '.')
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'is_hallucination': bool - True if predicted as hallucination
        - 'probability': float - Probability of hallucination (0-1)
        - 'method': str - Always 'ml_model'
        - 'error': Optional[str] - Error message if prediction failed
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    
    # Check if dependencies are available
    if not _import_dependencies():
        return {
            'is_hallucination': None,
            'probability': None,
            'method': 'ml_model',
            'error': 'Required dependencies (deepchem, xgboost) are not installed'
        }
    
    # Try to load model if not already loaded
    success, error_msg = load_ml_model()
    if not success:
        return {
            'is_hallucination': None,
            'probability': None,
            'method': 'ml_model',
            'error': error_msg
        }
    
    # Validate SMILES
    if not is_valid_smiles(product_smiles):
        return {
            'is_hallucination': True,  # Invalid SMILES treated as hallucination
            'probability': 1.0,
            'method': 'ml_model',
            'error': 'Invalid product SMILES'
        }
    
    if not is_valid_smiles(reactants_smiles):
        return {
            'is_hallucination': True,  # Invalid SMILES treated as hallucination
            'probability': 1.0,
            'method': 'ml_model',
            'error': 'Invalid reactant SMILES'
        }
    
    try:
        # Featurize product and reactants
        product_features = _featurizer.featurize([product_smiles])
        reactant_features = _featurizer.featurize([reactants_smiles])
        
        # Combine features (same as training)
        combined_features = np.concatenate([product_features, reactant_features], axis=1)
        
        # Get probability predictions
        probabilities = _ml_model.predict_proba(combined_features)[0]
        hallucination_prob = probabilities[1]  # Probability of class 1 (hallucination)
        
        # Use optimal threshold for binary prediction
        is_hallucination = hallucination_prob >= _optimal_threshold
        
        return {
            'is_hallucination': bool(is_hallucination),
            'probability': float(hallucination_prob),
            'method': 'ml_model',
            'error': None
        }
        
    except Exception as e:
        error_msg = f"Error during ML prediction: {str(e)}"
        log_message(f"ML Model: {error_msg}", logger)
        return {
            'is_hallucination': None,
            'probability': None,
            'method': 'ml_model',
            'error': error_msg
        }


def ml_hallucination_checker(product: str, res_smiles: list) -> Tuple[int, list]:
    """
    ML-based hallucination checker wrapper function.
    
    Similar interface to rule-based hallucination_checker() but uses ML model.
    Filters pathways based on ML predictions.
    
    Parameters
    ----------
    product : str
        SMILES string of the product molecule
    res_smiles : list
        List of list of reactant SMILES strings
    
    Returns
    -------
    Tuple[int, list]
        (status_code, valid_pathways)
        - status_code: 200 if successful, 500 if error
        - valid_pathways: List of pathways that passed ML hallucination check
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    valid_pathways = []
    
    for idx, smile_list in enumerate(res_smiles):
        if isinstance(smile_list, list):
            # Combine multiple reactants with '.'
            smiles_combined = ".".join(smile_list)
            
            if not is_valid_smiles(smiles_combined):
                log_message(f"ML Model: Invalid SMILES string: {smiles_combined}", logger)
                continue  # Skip invalid SMILES
            
            # Get ML prediction
            prediction = predict_hallucination_ml(product, smiles_combined)
            
            if prediction.get('error'):
                log_message(f"ML Model: Error predicting for pathway {idx}: {prediction['error']}", logger)
                # On error, be conservative and skip this pathway
                continue
            
            # If not hallucination, keep the pathway
            if not prediction['is_hallucination']:
                valid_pathways.append(smile_list)
                log_message(
                    f"ML Model: Pathway {idx} passed (probability={prediction['probability']:.3f})",
                    logger
                )
            else:
                log_message(
                    f"ML Model: Pathway {idx} filtered out (hallucination probability={prediction['probability']:.3f})",
                    logger
                )
        else:
            # Handle single SMILES string (shouldn't happen in normal flow, but handle gracefully)
            if is_valid_smiles(smile_list):
                prediction = predict_hallucination_ml(product, smile_list)
                
                if prediction.get('error'):
                    log_message(f"ML Model: Error predicting: {prediction['error']}", logger)
                    continue
                
                if not prediction['is_hallucination']:
                    valid_pathways.append([smile_list])
    
    log_message(f"ML Model: Valid pathways after filtering: {len(valid_pathways)}", logger)
    return 200, valid_pathways

