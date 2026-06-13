"""ML-based hallucination detection using trained XGBoost model."""

import os
import io
import json
import pickle
from typing import Dict, Optional, Tuple

import boto3
import numpy as np

from src.utils.ml_utils import extract_domain_features_single
from src.utils.utils_molecule import is_valid_smiles
from src.utils.job_context import logger as context_logger

_S3_BUCKET = os.getenv("ML_MODEL_S3_BUCKET", "recursivellm")
_S3_PREFIX = os.getenv("ML_MODEL_S3_PREFIX", "models")
ENABLE_LOGGING = os.getenv("ENABLE_LOGGING", "true").lower() != "false"


def _get_logger():
    """
    Retrieve the context logger when available.

    Returns
    -------
    logger or None
        Bound structlog logger, or None if logging is disabled or
        the context variable has not been set.
    """
    if not ENABLE_LOGGING:
        return None
    try:
        return context_logger.get()
    except LookupError:
        return None


def _log(message, logger=None):
    """
    Log a message via the context logger, falling back to print.

    Parameters
    ----------
    message : str
        Text to log.
    logger : optional
        Bound structlog logger. Uses ``print`` when None.
    """
    if logger:
        logger.info(message)
    else:
        print(message)


class HallucinationModel:
    """
    Lazy-loaded XGBoost hallucination detector.

    On first use, downloads model artifacts from S3 directly into memory
    (no local file caching). Kept in memory for subsequent calls within
    the same process.
    """

    def __init__(self):
        self._model = None
        self._featurizer = None
        self._threshold: float = 0.5
        self._use_domain: bool = True

    def _download_s3(self, filename):
        """
        Download a file from S3 into memory using AWS credentials.

        Parameters
        ----------
        filename : str
            Object key under ``_S3_PREFIX`` (e.g. 'xgboost_model.pkl').

        Returns
        -------
        data : bytes
            Raw file contents.
        """
        s3 = boto3.client(
            's3',
            aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
            aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        )
        buf = io.BytesIO()
        s3.download_fileobj(_S3_BUCKET, f"{_S3_PREFIX}/{filename}", buf)
        buf.seek(0)
        return buf

    def load(self):
        """
        Load model, featurizer, and metadata from S3 into memory.

        Returns
        -------
        success : bool
        error : str or None
            Human-readable error message on failure.
        """
        if self._model is not None:
            return True, None

        try:
            _log("ML Model: Downloading from S3...", _get_logger())

            self._model = pickle.load(self._download_s3("xgboost_hallucination_model.pkl"))
            self._featurizer = pickle.load(self._download_s3("xgboost_featurizer.pkl"))
            meta = json.load(self._download_s3("xgboost_metadata.json"))

            self._threshold = meta.get('optimal_threshold', 0.5)
            self._use_domain = meta.get('domain_features', True)
            _log(f"ML Model: Loaded (threshold={self._threshold:.4f})", _get_logger())
            return True, None
        except Exception as e:
            self._model = None
            return False, f"Error loading ML model: {e}"

    def predict(self, product_smiles, reactants_smiles):
        """
        Predict whether a single reaction is hallucinated.

        Parameters
        ----------
        product_smiles : str
            SMILES of the target product.
        reactants_smiles : str
            SMILES of the proposed reactants (dot-separated when multiple).

        Returns
        -------
        result : dict
            Keys: is_hallucination (bool | None), probability (float | None),
            method ('ml_model'), error (str | None).
        """
        def _error(msg):
            return {'is_hallucination': None, 'probability': None,
                    'method': 'ml_model', 'error': msg}

        success, err = self.load()
        if not success:
            return _error(err)

        if not is_valid_smiles(product_smiles):
            return {'is_hallucination': True, 'probability': 1.0,
                    'method': 'ml_model', 'error': 'Invalid product SMILES'}
        if not is_valid_smiles(reactants_smiles):
            return {'is_hallucination': True, 'probability': 1.0,
                    'method': 'ml_model', 'error': 'Invalid reactant SMILES'}

        try:
            product_fp = self._featurizer.featurize([product_smiles])
            reactant_fp = self._featurizer.featurize([reactants_smiles])
            features = np.concatenate([product_fp, reactant_fp], axis=1)

            if self._use_domain:
                domain = extract_domain_features_single(product_smiles, reactants_smiles)
                features = np.concatenate([features, domain.reshape(1, -1)], axis=1)

            prob = self._model.predict_proba(features)[0][1]
            return {
                'is_hallucination': bool(prob >= self._threshold),
                'probability': float(prob),
                'method': 'ml_model',
                'error': None,
            }
        except Exception as e:
            _log(f"ML Model: Prediction error: {e}", _get_logger())
            return _error(f"Prediction error: {e}")

    def check_pathways(self, product, res_smiles):
        """
        Filter a list of pathways, keeping only non-hallucinated ones.

        Fails open: if the model cannot load, all pathways are passed
        through unchanged.

        Parameters
        ----------
        product : str
            Target product SMILES.
        res_smiles : list
            Each element is a reactant SMILES string or a list of strings.

        Returns
        -------
        status_code : int
            200 on success.
        valid_pathways : list
            Pathways that the model did not flag as hallucinated.
        """
        logger = _get_logger()
        _log(f"=== ML Hallucination Check: {len(res_smiles)} pathways for {product} ===", logger)

        success, err = self.load()
        if not success:
            _log(f"ML Model: Unavailable ({err}), passing all pathways", logger)
            return 200, res_smiles

        valid_pathways = []
        for idx, smile_list in enumerate(res_smiles):
            reactants = smile_list if isinstance(smile_list, list) else [smile_list]
            smiles_combined = ".".join(reactants)

            prediction = self.predict(product, smiles_combined)

            if prediction.get('error') and prediction['is_hallucination'] is None:
                _log(f"  Pathway {idx}: SKIP (error: {prediction['error']})", logger)
                continue

            if not prediction['is_hallucination']:
                valid_pathways.append(reactants)
                _log(f"  Pathway {idx}: PASS (prob={prediction['probability']:.3f})", logger)
            else:
                _log(f"  Pathway {idx}: REJECT (prob={prediction['probability']:.3f})", logger)

        _log(f"=== ML Result: {len(valid_pathways)}/{len(res_smiles)} passed ===", logger)
        return 200, valid_pathways


# Module-level singleton
_model = HallucinationModel()

# Public API — unchanged signatures for existing callers
predict_hallucination_ml = _model.predict
ml_hallucination_checker = _model.check_pathways
