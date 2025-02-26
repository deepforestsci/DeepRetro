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
         llm: str = "claude-3-opus-20240229",
         az_model: str = "USPTO",
         stability_flag: str = "False",
         hallucination_check: str = "False") -> Any:
    """Run the retrosynthesis on specific molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    llm : str, optional
        LLM model, by default "claude-3-opus-20240229"

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

    res = run_prithvi(molecule=smiles,
                      llm=llm,
                      az_model=az_model,
                      stability_flag=stability_flag,
                      hallucination_check=hallucination_check)
    logging.info(f"Retrosynthesis result: {res}")
    return res
