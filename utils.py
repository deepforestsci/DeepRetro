"""Utils for Retrosynthesis"""
import ast
import os
from typing import Any, Dict, List, Optional, Sequence
import time
from rdkit import Chem
from anthropic import Anthropic
from aizynthfinder.aizynthfinder import AiZynthFinder
import logging
from variables import BASIC_MOLECULES, SYS_PROMPT, USER_PROMPT, ENCODING_SCALABILITY, REACTION_ENCODING_NAMES
from variables import bcolors
from utils_molecule import validity_check, calc_chemical_formula, calc_mol_wt
from utils_molecule import get_reaction_type
from cache import cache_results
from dotenv import load_dotenv

# load environment variables
load_dotenv()
ANTHROPIC_API_KEY = os.getenv('ANTHROPIC_API_KEY')
AZ_MODEL_CONFIG_PATH = os.getenv('AZ_MODEL_CONFIG_PATH')

# Set up logging
# log file based on time
# make a folder for each day
if not os.path.exists('logs'):
    os.makedirs('logs')

if not os.path.exists(f'logs/{time.strftime("%Y-%m-%d")}'):
    os.makedirs(f'logs/{time.strftime("%Y-%m-%d")}')

logging.basicConfig(
    filename=
    f"logs/{time.strftime('%Y-%m-%d')}/{time.strftime('%H-%M-%S')}.log",
    level=logging.INFO,  # Log level (INFO will log both inputs and outputs)
    format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
)
logger = logging.getLogger(__name__)
logger.info("Logging started")

client = Anthropic(api_key=ANTHROPIC_API_KEY)

# Check if the configuration file exists in `aizynthfinder/models/config.yml`, else check `../aizynthfinder/models/config.yml` if still not foubnd, raise an error
try:
    with open("aizynthfinder/models/config.yml", "r") as file:
        config_filename = "aizynthfinder/models/config.yml"
        logger.info(
            f"Configuration file found in `aizynthfinder/models/config.yml`")
except FileNotFoundError:
    with open("../aizynthfinder/models/config.yml", "r") as file:
        config_filename = "../aizynthfinder/models/config.yml"
        logger.info(
            f"Configuration file found in `../aizynthfinder/models/config.yml`"
        )
except FileNotFoundError:
    logging.error(
        f"{bcolors.FAIL}Configuration file not found in `aizynthfinder/models/config.yml` or `../aizynthfinder/models/config.yml`{bcolors.ENDC}"
    )
    raise FileNotFoundError(
        f"Configuration file not found in `aizynthfinder/models/config.yml` or `../aizynthfinder/models/config.yml`"
    )

# TODO: add try except block for invalid smiles in the entire pipeline


@cache_results
def call_LLM(molecule: str,
             LLM: str = "claude-3-opus-20240229",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None):
    """Calls the LLM model to predict the next step"""

    logger.info(f"Calling LLM with molecule: {molecule}")
    if messages is None:
        messages = [{
            "role":
            "user",
            "content": [{
                'type':
                "text",
                "text":
                USER_PROMPT.replace('{target_smiles}', molecule)
            }]
        }]
    try:
        message = client.messages.create(model=LLM,
                                         max_tokens=1024,
                                         temperature=temperature,
                                         system=SYS_PROMPT,
                                         messages=messages,
                                         top_p=0.9)
        res_text = message.content[0].text
    except Exception as e:
        logger.info(f"Error in calling LLM: {e}")
        message = client.messages.create(model=LLM,
                                         max_tokens=1024,
                                         temperature=temperature,
                                         system=SYS_PROMPT,
                                         messages=messages,
                                         top_p=0.9)
        res_text = message.content[0].text
    logger.info(f"Received response from LLM: {res_text}")

    try:
        result_list = ast.literal_eval(res_text)
    except Exception as e:
        logger.info(f"Error in parsing response: {e}")
        raise ValueError("Please Retry")
    res_molecules = result_list['data']
    res_explanations = result_list['explanation']
    res_confidence = result_list['confidence_scores']

    return res_molecules, res_explanations, res_confidence


def llm_pipeline(molecule: str,
                 LLM: str = "claude-3-opus-20240229",
                 messages: Optional[list[dict]] = None):
    """Pipeline to call LLM and validate the results"""
    output_pathways = []
    run = 0.0
    while (output_pathways == [] and run < 0.6):
        logger.info(
            f"Running LLM for the {run}th time for molecule: {molecule}")
        res_molecules, res_explanations, res_confidence = call_LLM(
            molecule, LLM, run, messages)
        output_pathways, output_explanations, output_confidence = validity_check(
            molecule, res_molecules, res_explanations, res_confidence)
        run += 0.2
    # mapped = []
    # for idx, pathway in enumerate(output_pathways):
    #     mapped.append(
    #         [pathway, output_explanations[idx], output_confidence[idx]])
    # # Sort the pathways by molecular weight. sort the explanations and confidence in the same order
    # output_pathways, output_explanations, output_confidence = zip(
    #     *sorted(zip(output_pathways, output_explanations, output_confidence),
    #             key=lambda x: calc_mol_wt(x[0][0])))

    output_pathways.sort(
        key=lambda x: Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(x[0])))
    # TODO: fix the matching of output pathways and confidence
    return output_pathways, output_explanations, output_confidence


def llm_pipeline_alt(molecule: str,
                     LLM: str = "claude-3-opus-20240229",
                     messages: Optional[list[dict]] = None):
    """Pipeline to call LLM and validate the results"""
    # output_pathways = []
    run = 0.0
    res_pathways, res_explanations = call_LLM(molecule, LLM, run / 10,
                                              messages)

    output_pathways = validity_check(molecule, res_pathways)
    run += 0.1

    return output_pathways, res_explanations


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
                    "name":
                    "",
                    "chemical_formula":
                    CalcMolFormula(Chem.MolFromSmiles(data['smiles'])),
                    "mass":
                    ExactMolWt(Chem.MolFromSmiles(data['smiles']))
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
    _, type = get_reaction_type(mol1, mol2)
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
