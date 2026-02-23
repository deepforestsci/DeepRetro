"""Utils for Retrosynthesis"""
import os
from typing import Any, Dict, List, Optional, Sequence
import time
from rdkit import Chem
from aizynthfinder.aizynthfinder import AiZynthFinder
import logging
import rootutils
import structlog
from dotenv import load_dotenv

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

from src.variables import BASIC_MOLECULES
from src.prithvi import run_prithvi
from src.utils.custom_logging import setup_logging

# load environment variables
load_dotenv()


def main(smiles: str,
         llm: str = "claude-opus-4-20250514",
         az_model: str = "USPTO",
         stability_flag: str = "False",
         hallucination_check: str = "False",
         use_protecting_group_feature: bool = False,
         protecting_group_mode: str = "auto",
         protecting_group_selections: dict = None,
         protection_state=None) -> Any:
    """Run the retrosynthesis on specific molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    llm : str, optional
        LLM model, by default "claude-opus-4-20250514"
    az_model : str, optional
        AZ model to use, by default "USPTO"
    stability_flag : str, optional
        Enable stability checking, by default "False"
    hallucination_check : str, optional
        Enable hallucination checking, by default "False"
    use_protecting_group_feature : bool, optional
        Enable protecting group feature, by default False
    protecting_group_mode : str, optional
        Mode for protecting group selection: "auto" or "hitl", by default "auto"
    protecting_group_selections : dict, optional
        User-selected protecting groups keyed by site_id (HITL mode), by default None
    protection_state : ProtectionState, optional
        Tracks active protecting groups across the synthesis tree.

    Returns
    -------
    Any
        Returns result of retrosynthesis.
    """
    # Initialize generic logging
    setup_logging()

    log = structlog.get_logger()
    log.info("-" * 50)
    log.info("Application initialization complete")

    res = run_prithvi(
        molecule=smiles,
        llm=llm,
        az_model=az_model,
        stability_flag=stability_flag,
        hallucination_check=hallucination_check,
        use_protecting_group_feature=use_protecting_group_feature,
        protecting_group_mode=protecting_group_mode,
        protecting_group_selections=protecting_group_selections,
        protection_state=protection_state)
    logging.info(f"Retrosynthesis result: {res}")
    return res
