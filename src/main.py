"""Utils for Retrosynthesis"""
import os
from typing import Any, Dict, List, Optional, Sequence
import time
from rdkit import Chem
from aizynthfinder.aizynthfinder import AiZynthFinder
import logging
import rootutils
from dotenv import load_dotenv

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

from src.variables import BASIC_MOLECULES

from src.utils.llm import llm_pipeline
from src.cache import cache_results
from src.utils.parse import format_output

# load environment variables
load_dotenv()

# setup rootutils

AZ_MODEL_CONFIG_PATH = f"{root_dir}/{os.getenv('AZ_MODEL_CONFIG_PATH')}"

# Set up logging
# Create logs directory if it doesn't exist
if not os.path.exists(f'{root_dir}/logs'):
    os.makedirs(f'{root_dir}/logs')

# Create date-based logs directory if it doesn't exist
date_dir = f'{root_dir}/logs/{time.strftime("%Y-%m-%d")}'
if not os.path.exists(date_dir):
    os.makedirs(date_dir)

# Configure root logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)

try:
    with open(AZ_MODEL_CONFIG_PATH, "r") as file:
        logger.info(f"AZ_MODEL_CONFIG_PATH found: {AZ_MODEL_CONFIG_PATH}")
        config_filename = AZ_MODEL_CONFIG_PATH
except FileNotFoundError:
    logger.error(f"AZ_MODEL_CONFIG_PATH not found at {AZ_MODEL_CONFIG_PATH}")
    raise FileNotFoundError(
        f"AZ_MODEL_CONFIG_PATH not found at {AZ_MODEL_CONFIG_PATH}")

# TODO: add try except block for invalid smiles in the entire pipeline


@cache_results
def run_az(smiles: str) -> tuple[Any, Sequence[Dict[str, Any]]]:
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
    if smiles in BASIC_MOLECULES:
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
    # images = finder.routes.images
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
    if smiles in BASIC_MOLECULES:
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
    images = finder.routes.images
    return status, result_dict, images


def calc_yield(mol1, mol2):
    """Calculate the yield of a reaction"""
    return "#"


# recursive function to run pritvi
def rec_run_prithvi(molecule, job_id):
    solved, result_dict = run_az(molecule)
    result_dict = result_dict[0]
    if not solved:

        logger.info(f"AZ failed for {molecule}, running LLM")
        out_pathways, out_explained, out_confidence = llm_pipeline(molecule)
        result_dict = {
            'type':
            'mol',
            'smiles':
            molecule,
            # 'confidence': out_confidence,
            "is_chemical":
            True,
            "in_stock":
            False,
            'children': [{
                "type": "reaction",
                "is_reaction": True,
                "metadata": {
                    "policy_probability": out_confidence,
                },
                "children": []
            }]
        }
        logger.info(f"LLM returned {out_pathways}")
        logger.info(f"LLM explained {out_explained}")
        for pathway in out_pathways:
            if isinstance(pathway, list):
                temp_stat = []
                for mol in pathway:
                    res, stat = rec_run_prithvi(mol, job_id)
                    if stat:
                        temp_stat.append(True)
                        result_dict['children'][0]['children'].append(res)
                logger.info(f"temp_stat: {temp_stat}")
                if all(temp_stat):
                    solved = True
            else:
                res, solved = rec_run_prithvi(pathway, job_id)
                result_dict['children'][0]['children'].append(res)
            if solved:
                logger.info('breaking')
                break
    else:
        logger.info(f"AZ solved {molecule}")
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved


def run_prithvi(molecule):
    # Generate a unique job ID using timestamp and a random suffix
    job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"

    # Set up file handler for this specific job
    file_handler = logging.FileHandler(f"{date_dir}/job_{job_id}.log")
    file_handler.setFormatter(
        logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)

    logger.info("Starting new synthesis job")

    try:
        result_dict, solved = rec_run_prithvi(molecule, job_id)
        output_data = format_output(result_dict)
        return output_data
    finally:
        # Clean up handlers
        logger.removeHandler(file_handler)
