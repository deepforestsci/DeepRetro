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
         use_agentic_pipeline: bool = True,
         max_tool_iterations: int = 5) -> Any:
    """Run the retrosynthesis on specific molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    llm : str, optional
        LLM model, by default "claude-opus-4-20250514"
    az_model : str, optional
        AZ model to be used, by default "USPTO"
    stability_flag : str, optional
        Stability flag, by default "False"
    hallucination_check : str, optional
        Hallucination check, by default "False"
    use_protecting_group_feature : bool, optional
        Whether to use protecting group feature, by default False
    use_agentic_pipeline : bool, optional
        Whether to use agentic pipeline with tools, by default True
    max_tool_iterations : int, optional
        Maximum number of tool iterations for agentic pipeline, by default 5

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
        use_agentic_pipeline=use_agentic_pipeline,
        max_tool_iterations=max_tool_iterations)
    logging.info(f"Retrosynthesis result: {res}")
    return res
