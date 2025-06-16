"""Main entry point for orchestrating the retrosynthesis process.

This module provides the `main()` function, which serves as the primary
coordinator for running a retrosynthesis task on a given molecule. It handles
application initialization steps like setting up logging and then invokes
the core processing logic, which is currently delegated to the `run_prithvi`
function from the `src.prithvi` module.

Environment variables are loaded using `dotenv` at the module level.
"""
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
    """Runs the end-to-end retrosynthesis process for a target molecule.

    This function initializes logging, and then calls `run_prithvi` to perform
    the retrosynthesis using a specified Large Language Model (LLM) and
    AiZynthFinder (AZ) model configuration. It can also enable optional
    stability and hallucination checks based on the provided flags.

    Args:
        smiles (str):
            The SMILES string of the target molecule for retrosynthesis.
        llm (str, optional):
            Identifier for the Large Language Model to be used. Defaults to
            "claude-3-opus-20240229". This is passed to `run_prithvi`.
        az_model (str, optional):
            Identifier for the AiZynthFinder model/configuration to use.
            Defaults to "USPTO". This is passed to `run_prithvi`.
        stability_flag (str, optional):
            Flag to enable/disable stability checks, passed to `run_prithvi`.
            Expected values are "True" or "False" (as strings). Defaults to "False".
        hallucination_check (str, optional):
            Flag to enable/disable hallucination checks, passed to `run_prithvi`.
            Expected values are "True" or "False" (as strings). Defaults to "False".

    Returns:
        Any:
            The result of the retrosynthesis process, as returned by `run_prithvi`.
            This is typically a dictionary containing the proposed synthetic routes,
            explanations, scores, and other metadata.
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
