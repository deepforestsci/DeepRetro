"""Metrics used by parse/format_output (ported from src.utils.utils_molecule subset)."""

from __future__ import annotations

import os
from pathlib import Path

import joblib

from deepretro.migration import job_context
from deepretro.utils import utils_molecule
from deepretro.utils.variables import ENCODING_SCALABILITY, REACTION_ENCODING_NAMES

ENABLE_LOGGING = os.getenv("ENABLE_LOGGING", "true").lower() != "false"


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    for d in [here, *here.parents]:
        if (d / "config" / "langfuse_config.json").is_file():
            return d
    raise FileNotFoundError(
        f"Could not locate repo root from {here} (expected config/langfuse_config.json)"
    )


def _rxn_classification_model_path() -> str:
    rel = os.getenv("RXN_CLASSIFICATION_MODEL_PATH", "") or ""
    if not rel.strip():
        return ""
    p = Path(rel)
    if p.is_absolute():
        return str(p)
    return str(_repo_root() / rel)


def log_message(message: str, logger=None) -> None:
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def calc_confidence_estimate(probability: float) -> float:
    """Confidence estimate from policy probability (same rules as src)."""
    if isinstance(probability, list):
        probability = probability[0] if probability else 0.0
    if probability < 0.3:
        probability = 1 - probability
    elif probability < 0.45 and probability >= 0.3:
        probability += 0.5
    elif probability < 0.6 and probability >= 0.45:
        probability += 0.3

    probability = round(probability, 2)
    if probability > 0.99:
        probability = 0.99
    return probability


def get_reaction_type(mol1, mol2, model_path: str):
    """Load sklearn joblib classifier and predict reaction encoding index."""
    try:
        clf = joblib.load(model_path)
        mol1_fingerprint = utils_molecule.compute_fingerprint(mol1)
        mol2_fingerprint = utils_molecule.compute_fingerprint(mol2)
        if mol1_fingerprint is None or mol2_fingerprint is None:
            return "Unknown Reaction", -1
        reaction_type = clf.predict([mol1_fingerprint + mol2_fingerprint])
        idx = int(reaction_type[0])
        return REACTION_ENCODING_NAMES[idx], idx
    except FileNotFoundError:
        return "Unknown Reaction", -1
    except Exception as e:
        logger = job_context.logger.get() if ENABLE_LOGGING else None
        log_message(f"Error in get_reaction_type: {e}", logger)
        return "Unknown Reaction", -1


def calc_scalability_index(mol1, mol2):
    """Scalability label from reaction type model; ``N/A`` if model missing or error."""
    path = _rxn_classification_model_path()
    if not path:
        return "N/A"
    try:
        _, rtype = get_reaction_type(mol1, mol2, path)
        if rtype == -1:
            return "N/A"
        return str(ENCODING_SCALABILITY[rtype])
    except Exception as e:
        logger = job_context.logger.get() if ENABLE_LOGGING else None
        log_message(f"Error in calc_scalability_index: {e}", logger)
        return "N/A"
