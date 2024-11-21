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


def main(smiles: str) -> Any:
    # Initialize generic logging
    setup_logging()

    log = structlog.get_logger()
    log.info("Application initialization complete")

    res = run_prithvi(smiles)
    logging.info(f"Retrosynthesis result: {res}")
    return res
