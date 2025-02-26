"""Utils for AiZynthFinder"""
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
    """Run the retrosynthesis using AiZynthFinder

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule

    Returns
    -------
    tuple[Any, Sequence[Dict[str, Any]]]
        A tuple containing the status of the retrosynthesis, 
        the results dictionary
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
    """Run the retrosynthesis using AiZynthFinder

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule

    Returns
    -------
    tuple[Any, Sequence[Dict[str, Any]]]
        A tuple containing the status of the retrosynthesis, 
        the results dictionary
        
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
    """Check if the molecule is a basic molecule
    (if number of C atons is less than 5)

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule

    Returns
    -------
    bool
        True if the molecule is a basic molecule, False otherwise
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
