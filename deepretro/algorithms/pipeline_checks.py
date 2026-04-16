"""LLM pathway filters: same behavior as src stability/hallucination wrappers, using deepretro algorithms."""

from __future__ import annotations

import os

from deepretro.algorithms.hallucination_checker import calculate_hallucination_score
from deepretro.algorithms.stability_checker import check_molecule_stability
from deepretro import logging as job_context
from deepretro.utils.utils_molecule import is_valid_smiles

ENABLE_LOGGING = os.getenv("ENABLE_LOGGING", "true").lower() != "false"


def log_message(message: str, logger=None) -> None:
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def stability_checker(res_smiles: list):
    """Filter pathways by heuristic stability (RDKit), matching src.utils.stability_checks API.

    Returns
    -------
    tuple[int, list]
        Status code 200 and list of pathways that pass stability thresholds.
    """
    logger = job_context.logger.get() if ENABLE_LOGGING else None
    valid_pathways = []
    for idx, smile_list in enumerate(res_smiles):
        valid = []
        if isinstance(smile_list, list):
            for smiles in smile_list:
                if not is_valid_smiles(smiles):
                    log_message("Invalid SMILES string", logger)

                stability_dict = check_molecule_stability(smiles)
                log_message(f"Stability dict: {stability_dict}", logger)

                if stability_dict["assessment"] in ("Likely stable", "Moderately stable"):
                    valid.append(smiles)
            if len(valid) == len(smile_list):
                valid_pathways.append(valid)
        else:
            if is_valid_smiles(smile_list):
                stability_dict = check_molecule_stability(smile_list)
                log_message(f"Stability dict: {stability_dict}", logger)

                if stability_dict["assessment"] in ("Likely stable", "Moderately stable"):
                    valid_pathways.append([smile_list])
    log_message(f"Valid pathways: {valid_pathways} after stability check", logger)
    return 200, valid_pathways


def hallucination_checker(product: str, res_smiles: list):
    """Filter pathways by hallucination score, matching src.utils.hallucination_checks API.

    Notes
    -----
    The single-SMILES branch passes ``(reactant, product)`` into
    ``calculate_hallucination_score``; src called with one argument only, which
    does not match the deepretro signature.
    """
    logger = job_context.logger.get() if ENABLE_LOGGING else None
    valid_pathways = []
    for idx, smile_list in enumerate(res_smiles):
        if isinstance(smile_list, list):
            smiles_combined = ".".join(smile_list)
            if not is_valid_smiles(smiles_combined):
                log_message(f"Invalid SMILES string: {smiles_combined}", logger)

            hallucination_report = calculate_hallucination_score(
                smiles_combined, product
            )
            log_message(f"Hallucination report: {hallucination_report}", logger)

            if hallucination_report["severity"] in ("low", "medium"):
                valid_pathways.append(smile_list)
        else:
            if is_valid_smiles(smile_list):
                hallucination_report = calculate_hallucination_score(
                    smile_list, product
                )
                log_message(f"Hallucination report: {hallucination_report}", logger)

                if hallucination_report["severity"] in ("low", "medium"):
                    valid_pathways.append([smile_list])
    log_message(f"Valid pathways: {valid_pathways}", logger)
    return 200, valid_pathways
