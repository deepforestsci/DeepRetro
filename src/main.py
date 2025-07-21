"""Utils for Retrosynthesis

This module provides the main entry point for running retrosynthesis analysis using DeepRetro.
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
         llm: str = "claude-opus-4-20250514",
         az_model: str = "USPTO",
         stability_flag: str = "False",
         hallucination_check: str = "False") -> Any:
    """Run the retrosynthesis on a specific molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    llm : str, optional
        LLM model to use, by default "claude-opus-4-20250514".
    az_model : str, optional
        Retrosynthesis model backend, by default "USPTO".
    stability_flag : str, optional
        Enable stability validation, by default "False".
    hallucination_check : str, optional
        Enable hallucination detection, by default "False".

    Returns
    -------
    Any
        Result of retrosynthesis.

    Examples
    --------
    >>> result = main("COc1ccc(-c2ccc(/C=C(\\C#N)c3ccc(-c4ccncc4)cc3)cc2)cc1")  # doctest: +SKIP
    >>> isinstance(result, dict)  # doctest: +SKIP
    True
    >>> # Example with custom model and validation
    >>> result = main(
    ...     "COc1ccc(-c2ccc(/C=C(\\C#N)c3ccc(-c4ccncc4)cc3)cc2)cc1",
    ...     llm="claude-3-sonnet",
    ...     az_model="Pistachio_100",
    ...     stability_flag="True",
    ...     hallucination_check="True"
    ... )  # doctest: +SKIP
    >>> isinstance(result, dict)  # doctest: +SKIP
    True
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
