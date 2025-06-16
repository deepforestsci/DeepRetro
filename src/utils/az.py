"""Utilities for interacting with AiZynthFinder for retrosynthesis tasks.

This module provides functions to run retrosynthetic analysis using AiZynthFinder,
including loading appropriate models and configurations. It also includes helper
functions, such as checking if a molecule is considered "basic" to potentially
skip expensive computations.

Results from AiZynthFinder runs are cached using the `cache_results` decorator
to improve performance on repeated requests for the same molecule.

Configuration for AiZynthFinder model paths (e.g., `AZ_MODEL_CONFIG_PATH`,
`AZ_MODELS_PATH`) is expected to be set via environment variables and is
resolved relative to the project root determined by `rootutils`.
"""
import os
from aizynthfinder.aizynthfinder import AiZynthFinder
from typing import Any, Dict, List, Optional, Sequence
from src.variables import BASIC_MOLECULES, ENCODING_SCALABILITY
from src.cache import cache_results
from src.utils.job_context import logger as context_logger
import rootutils
from rdkit import Chem
from rdkit.Chem import rdqueries

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

AZ_MODEL_CONFIG_PATH = f"{root_dir}/{os.getenv('AZ_MODEL_CONFIG_PATH')}"
AZ_MODELS_PATH = f"{root_dir}/{os.getenv('AZ_MODELS_PATH')}"


@cache_results
def run_az(smiles: str,
           az_model: str = "USPTO") -> tuple[Any, Sequence[Dict[str, Any]]]:
    """Run a retrosynthesis prediction using AiZynthFinder for a given SMILES.

    This function initializes AiZynthFinder with a specified configuration,
    attempts to find retrosynthetic routes for the target molecule, and returns
    the status and results.

    The function first tries to load a model-specific configuration file from
    `AZ_MODELS_PATH/<az_model>/config.yml`. If not found, it falls back to
    the general `AZ_MODEL_CONFIG_PATH`.

    If the input SMILES string corresponds to a "basic" molecule (as determined
    by `is_basic_molecule()` or a predefined list `BASIC_MOLECULES`), the
    AiZynthFinder analysis is skipped, and a predefined "in_stock" result is
    returned.

    Args:
        smiles (str): The SMILES string of the target molecule for retrosynthesis.
        az_model (str, optional): The name of the AiZynthFinder model to use.
            This typically corresponds to a subdirectory within `AZ_MODELS_PATH`
            containing the model files and a `config.yml`. Defaults to "USPTO".

    Returns:
        tuple[bool, Sequence[Dict[str, Any]]]:
            A tuple containing:
            - status (bool): True if the retrosynthesis was successful (solved),
              False otherwise.
            - result_dict (Sequence[Dict[str, Any]]): A list of dictionaries
              representing the found routes, in the format provided by
              AiZynthFinder's `routes.dict_with_extra()` method. If the molecule
              is basic, this will be a predefined list indicating the molecule
              is in stock.

    Raises:
        FileNotFoundError: If the AiZynthFinder configuration file cannot be found
            after checking both model-specific and general paths.
    """
    try:
        config_path = f"{AZ_MODELS_PATH}/{az_model}/config.yml"
        with open(config_path, "r") as file:
            logger = context_logger.get()
            logger.info(f"AZ_MODEL_CONFIG_PATH found: {config_path}")
            config_filename = config_path
    except FileNotFoundError:
        logger.error(f"AZ_MODEL_CONFIG_PATH not found at {config_path}")
        try:
            with open(AZ_MODEL_CONFIG_PATH, "r") as file:
                logger = context_logger.get()
                logger.info(
                    f"AZ_MODEL_CONFIG_PATH found: {AZ_MODEL_CONFIG_PATH}")
                config_filename = AZ_MODEL_CONFIG_PATH
        except FileNotFoundError:
            logger.error(
                f"AZ_MODEL_CONFIG_PATH not found at {AZ_MODEL_CONFIG_PATH}")
            raise FileNotFoundError(
                f"AZ_MODEL_CONFIG_PATH not found at {AZ_MODEL_CONFIG_PATH}")
    # if simple molecule, skip the retrosynthesis
    if smiles in BASIC_MOLECULES or is_basic_molecule(smiles):
        return True, [{
            'type': 'mol',
            'hide': False,
            'smiles': smiles,
            'is_chemical': True,
            'in_stock': True,
        }]
    finder = AiZynthFinder(configfile=config_filename)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats = finder.extract_statistics()
    status = stats['is_solved']
    result_dict = finder.routes.dict_with_extra(include_metadata=True,
                                                include_scores=True)
    return status, result_dict


@cache_results
def run_az_with_img(smiles: str) -> tuple[Any, Sequence[Dict[str, Any]]]:
    """Runs retrosynthesis using AiZynthFinder and extracts route images.

    Similar to `run_az`, this function performs retrosynthetic analysis but also
    extracts graphical representations (images) of the resulting routes.
    It uses a fixed configuration path defined by `AZ_MODEL_CONFIG_PATH`.

    If the input SMILES string corresponds to a "basic" molecule (as determined
    by `is_basic_molecule()` or a predefined list `BASIC_MOLECULES`), the
    AiZynthFinder analysis is skipped, and a predefined "in_stock" result is
    returned (with no images).

    Args:
        smiles (str): The SMILES string of the target molecule.

    Returns:
        tuple[bool, Sequence[Dict[str, Any]], Optional[List[Any]]]:
            A tuple containing:
            - status (bool): True if retrosynthesis was successful (solved),
              False otherwise.
            - result_dict (Sequence[Dict[str, Any]]): A list of dictionaries
              representing the found routes, from `routes.dict_with_extra()`.
              If basic, a predefined "in_stock" list.
            - images (Optional[List[Any]]): A list of images representing the
              routes, as provided by `finder.routes.images`. This will be `None`
              or an empty list if the molecule is basic or no routes are found.
              The exact type of items in the list depends on AiZynthFinder's
              image generation (e.g., file paths, image data).

    Raises:
        FileNotFoundError: If the `AZ_MODEL_CONFIG_PATH` is not found.
    """
    # if simple molecule, skip the retrosynthesis
    if smiles in BASIC_MOLECULES or is_basic_molecule(smiles):
        return True, [{
            'type': 'mol',
            'hide': False,
            'smiles': smiles,
            'is_chemical': True,
            'in_stock': True,
        }]
    finder = AiZynthFinder(configfile=AZ_MODEL_CONFIG_PATH)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats = finder.extract_statistics()
    status = stats['is_solved']
    result_dict = finder.routes.dict_with_extra(include_metadata=True,
                                                include_scores=True)
    images = finder.routes.images
    return status, result_dict, images


def is_basic_molecule(smiles: str) -> bool:
    """Checks if a molecule is considered "basic" for retrosynthesis purposes.

    A molecule is considered basic if it meets either of the following criteria,
    determined using RDKit:
    1. The total number of atoms in the molecule is less than 5.
    2. The number of carbon atoms in the molecule is less than 5.

    This check is used to quickly identify simple molecules for which running
    a full retrosynthesis analysis might be unnecessary, as they might be
    considered readily available or too simple to break down further.

    Args:
        smiles (str): The SMILES string of the molecule to check.

    Returns:
        bool:
            True if the molecule is considered basic, False otherwise.
            Returns False if the SMILES string cannot be parsed by RDKit.
    """
    #
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False

    q = rdqueries.AtomNumEqualsQueryAtom(6)
    num_c_atoms = len(mol.GetAtomsMatchingQuery(q))
    # if total number of atoms is less than 5, return True
    if mol.GetNumAtoms() < 5:
        return True
    elif num_c_atoms < 5:
        return True
    # if total number of C atoms is less than 5, return True

    return False
