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

from src.variables import BASIC_MOLECULES, ENCODING_SCALABILITY
from src.variables import bcolors
from src.utils.utils_molecule import validity_check, calc_chemical_formula, calc_mol_wt
from src.utils.utils_molecule import get_reaction_type
from src.utils.llm import call_LLM, llm_pipeline
from src.cache import cache_results

# load environment variables
load_dotenv()

# setup rootutils

AZ_MODEL_CONFIG_PATH = f"{root_dir}/{os.getenv('AZ_MODEL_CONFIG_PATH')}"
RXN_CLASSIFICATION_MODEL_PATH = f"{root_dir}/{os.getenv('RXN_CLASSIFICATION_MODEL_PATH')}"

# Set up logging
# log file based on time
# make a folder for each day
if not os.path.exists(f'{root_dir}/logs'):
    os.makedirs(f'{root_dir}/logs')

if not os.path.exists(f'{root_dir}/logs/{time.strftime("%Y-%m-%d")}'):
    os.makedirs(f'{root_dir}/logs/{time.strftime("%Y-%m-%d")}')

logging.basicConfig(
    filename=
    f"{root_dir}/logs/{time.strftime('%Y-%m-%d')}/{time.strftime('%H-%M-%S')}.log",
    level=logging.INFO,  # Log level (INFO will log both inputs and outputs)
    format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
)
logger = logging.getLogger(__name__)
logger.info("Logging started")

# Check if the configuration file exists in `aizynthfinder/models/config.yml`, else check `../aizynthfinder/models/config.yml` if still not foubnd, raise an error
try:
    with open(f"{root_dir}/aizynthfinder/models/config.yml", "r") as file:
        config_filename = f"{root_dir}/aizynthfinder/models/config.yml"
        logger.info(
            f"Configuration file found in `aizynthfinder/models/config.yml`")
except FileNotFoundError:
    try:
        with open(f"{root_dir}/../aizynthfinder/models/config.yml",
                  "r") as file:
            config_filename = f"{root_dir}/../aizynthfinder/models/config.yml"
            logger.info(
                f"Configuration file found in `../aizynthfinder/models/config.yml`"
            )
    except FileNotFoundError:
        logging.error(
            f"Configuration file not found in `aizynthfinder/models/config.yml` or `../aizynthfinder/models/config.yml`"
        )
        raise FileNotFoundError(
            f"Configuration file not found in `aizynthfinder/models/config.yml` or `../aizynthfinder/models/config.yml`"
        )

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


def parse_step(data,
               include_metadata: Optional[bool] = None,
               step_list: Optional[List] = None,
               dependency_list: Optional[Dict] = None,
               parent_id: Optional[int] = None) -> Dict:
    """
    Parses the input data recursively to extract steps and their dependencies in a specified format.

    Parameters
    ----------
    data: Dict
        A dictionary containing the SMILES representation and optional children of a molecule.
    include_metadata: bool, Optional
        A flag to include metadata in the output. Defaults to None, which includes metadata.
    step_list: List, Optional
        A list of steps extracted so far. Defaults to None, which initializes an empty list.
    dependency_list: Dict, Optional
        A dictionary mapping parent step IDs to their child step IDs. Defaults to None, which 
        initializes an empty dictionary.
    parent_id: int, Optional
        The ID of the parent step. Defaults to None for the root node.

    Returns
    -------
    Dict
        A dictionary containing:
            - 'dependencies': Dict
                A mapping of step IDs to their respective child step IDs.
            - 'steps': List[Dict]
                A list of steps where each step is represented by a dictionary containing:
                    - 'step': str
                        The identifier for the step.
                    - 'reactants': List
                        A list of reactants for the step, each represented by:
                            - 'smiles': str
                                SMILES representation of the reactant.
                            - 'reactant_metadata': Any
                                Metadata related to the reactant (currently not populated).
                    - 'reagents': List
                        A list of reagents for the step (currently empty).
                    - 'products': List
                        A list of products for the step, each represented by:
                            - 'smiles': str
                                SMILES representation of the product.
                            - 'product_metadata': Any
                                Metadata related to the product (currently not populated).
                    - 'conditions': List
                        A list of conditions for the step (currently empty).
                    - 'reactionmetrics': List
                        A list of reaction metrics for the step (currently empty).
    """
    if step_list is None:
        step_list = []
    if dependency_list is None:
        dependency_list = {}
    if include_metadata is None:
        include_metadata = True

    step_id = len(step_list) + 1

    if 'children' in data:
        step = {
            "step":
            str(step_id),
            "reactants": [],
            "reagents": [],
            "products": [{
                "smiles": data['smiles'],
                "product_metadata": {
                    "name": "",
                    "chemical_formula": calc_chemical_formula(data['smiles']),
                    "mass": calc_mol_wt(data['smiles'])
                }
            }],
            "conditions": [],
            "reactionmetrics": [{
                "scalabilityindex":
                "",
                "confidenceestimate":
                calc_confidence_estimate(
                    data['children'][0]['metadata']['policy_probability']),
                "closestliterature":
                "",
                # "yield":
                # ""
            }]
        }
        if parent_id is not None:

            if str(parent_id) in dependency_list:
                # print(f"Parent ID is dependency list: {parent_id}")
                dependency_list[str(parent_id)].append(str(step_id))
            else:
                # print(f"Parent ID is not in dependency list: {parent_id}")
                dependency_list[str(parent_id)] = [str(step_id)]
    else:

        step = None
        dependency_list[str(parent_id)] = []
    if parent_id is not None and not data.get('is_reaction', False):
        if data['smiles'] in BASIC_MOLECULES:
            step_list[parent_id - 1]['reagents'].append({
                "smiles":
                data['smiles'],
                "reagent_metadata": {
                    "name": "",
                    "chemical_formula": calc_chemical_formula(data['smiles']),
                    "mass": calc_mol_wt(data['smiles'])
                }
            })
        else:
            step_list[parent_id - 1]["reactants"].append({
                "smiles":
                data['smiles'],
                "reactant_metadata": {
                    "name": "",
                    "chemical_formula": calc_chemical_formula(data['smiles']),
                    "mass": calc_mol_wt(data['smiles'])
                }
            })
        # step_list[parent_id - 1]["reactionmetrics"][0]["yield"] = calc_yield(
        #     data['smiles'], step_list[parent_id - 1]["products"][0]["smiles"])

        step_list[parent_id - 1]["reactionmetrics"][0][
            "scalabilityindex"] = calc_scalability_index(
                data['smiles'],
                step_list[parent_id - 1]["products"][0]["smiles"])
        # step_list[parent_id -
        #           1]["reactionmetrics"][0]["reactiontype"] = get_reaction_type(
        #               data['smiles'],
        #               step_list[parent_id - 1]["products"][0]["smiles"])[0]

    if step is not None:
        step_list.append(step)
        if 'children' in data:
            for child in data['children']:
                if 'children' in child:
                    for c in child['children']:
                        parse_step(c, include_metadata, step_list,
                                   dependency_list, step_id)
                    # parse_step(child['children'], step_list, dependency_list, step_id)
                else:
                    parse_step(child, include_metadata, step_list,
                               dependency_list, step_id)
    # print(f"Step ID: {step_id}")
    # print(f"Parent ID: {parent_id}")
    # print(f"dependency_list: {dependency_list}")
    # print(f"Data: {data}")
    # print("--------------------")
    return {"dependencies": dependency_list, "steps": step_list}


def format_output(data):
    """Format the output data for visualization"""
    output_data = parse_step(data)
    # fix_dependencies(output_data['dependencies'], output_data['steps'])
    output_data['dependencies'] = fix_dependencies(output_data['dependencies'],
                                                   output_data['steps'])
    return output_data


def fix_dependencies(dependencies, step_list):
    """Fix the dependencies to be in the correct format"""
    fixed_dependencies = {}
    storage = {}
    for dicti in step_list:
        storage[dicti['products'][0]['smiles']] = dicti['step']
        fixed_dependencies[dicti['step']] = []
    for dicti in step_list:
        for react in dicti['reactants']:
            if react['smiles'] in storage:
                if dicti['step'] in fixed_dependencies:
                    fixed_dependencies[dicti['step']].append(
                        storage[react['smiles']])
                else:
                    fixed_dependencies[dicti['step']] = [
                        storage[react['smiles']]
                    ]

        pass
    return fixed_dependencies


def calc_yield(mol1, mol2):
    """Calculate the yield of a reaction"""
    return "#"


def calc_scalability_index(mol1, mol2):
    """Calculate the scalability index of a reaction"""
    _, type = get_reaction_type(mol1, mol2, RXN_CLASSIFICATION_MODEL_PATH)
    return str(ENCODING_SCALABILITY[type])


def calc_confidence_estimate(probability: float) -> float:
    """Calculate the confidence estimate based on the probability

    Parameters
    ----------
    probability : float
        Probability of the prediction

    Returns
    -------
    float
        Confidence estimate
    """
    if isinstance(probability, list):
        probability = probability[0]
    if probability < 0.3:
        probability = 1 - probability
    elif probability < 0.45 and probability >= 0.3:
        probability += 0.5
    elif probability < 0.6 and probability >= 0.45:
        probability += 0.3

    # limit the confidence estimate to 2 decimal places, round to the
    # nearest 0.01
    probability = round(probability, 2)
    if probability > 0.99:
        probability = 0.99
    return probability


# recursive function to run pritvi
def rec_run_prithvi(molecule):
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
                    res, stat = rec_run_prithvi(mol)
                    if stat:
                        temp_stat.append(True)
                        result_dict['children'][0]['children'].append(res)
                logger.info(f"temp_stat: {temp_stat}")
                if all(temp_stat):
                    solved = True
            else:
                res, solved = rec_run_prithvi(pathway)
                result_dict['children'][0]['children'].append(res)
            if solved:
                logger.info('breaking')
                break
    else:
        logger.info(f"AZ solved {molecule}")
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved


def run_prithvi(molecule):
    result_dict, solved = rec_run_prithvi(molecule)
    output_data = format_output(result_dict)
    return output_data
