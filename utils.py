"""Utils for Retrosynthesis"""
import ast
from typing import Any, Dict, List, Optional, Sequence
import time
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from anthropic import Anthropic
from aizynthfinder.aizynthfinder import AiZynthFinder
import diskcache as dc
import logging
from variables import BASIC_MOLECULES, SYS_PROMPT, USER_PROMPT, ENCODING_SCALABILITY, REACTION_ENCODING_NAMES
from variables import bcolors
import joblib

# Set up logging
# log file based on time
logging.basicConfig(
    filename=f"logs/{time.strftime('%Y-%m-%d-%H-%M-%S')}.log",
    level=logging.INFO,  # Log level (INFO will log both inputs and outputs)
    format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
)

# Create a disk cache
cache = dc.Cache('cache_dir_chall_failed')

client = Anthropic(
    api_key=
    "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",
)

config_filename = "../aizynthfinder/models/config.yml"


def clear_cache():
    """Clear the cache"""
    cache.clear()


def cache_results(func):
    """Decorator to cache results using diskcache"""

    def wrapper(*args, **kwargs):
        cache_key = func.__name__ + "_" + str(args) + str(
            kwargs)  # Unique key with function name
        if cache_key in cache:
            return cache[cache_key]
        else:
            result = func(*args, **kwargs)
            cache[cache_key] = result
            return result

    return wrapper


# TODO: add try except block for invalid smiles in the entire pipeline


def is_valid_smiles(smiles: str) -> bool:
    """Check if the SMILES string is valid

    Parameters
    ----------
    smiles : str
        smiles string

    Returns
    -------
    bool
        True if the smiles is valid, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False
    if mol is None:
        return False
    # Additional checks if needed
    return True


def substructure_matching(target_smiles: str, query_smiles: str):
    # Convert SMILES to RDKit molecule objects
    try:
        target_molecule = Chem.MolFromSmiles(target_smiles)
    except:
        logging.info(f"Error in parsing target molecule: {target_smiles}")

    try:
        query_molecule = Chem.MolFromSmiles(query_smiles)
    except:
        logging.info(f"Error in parsing query molecule: {query_smiles}")

    # Check if the query substructure is present in the target molecule
    try:
        if target_molecule.HasSubstructMatch(query_molecule):
            return 1
        else:
            return 0
    except:
        return 0


@cache_results
def call_LLM(molecule: str,
             LLM: str = "claude-3-opus-20240229",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None):
    """Calls the LLM model to predict the next step"""

    logging.info(f"Calling LLM with molecule: {molecule}")
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
        logging.info(f"Error in calling LLM: {e}")
        message = client.messages.create(model=LLM,
                                         max_tokens=1024,
                                         temperature=temperature,
                                         system=SYS_PROMPT,
                                         messages=messages,
                                         top_p=0.9)
        res_text = message.content[0].text
    logging.info(f"Received response from LLM: {res_text}")

    try:
        result_list = ast.literal_eval(res_text)
    except Exception as e:
        logging.info(f"Error in parsing response: {e}")
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
        logging.info(
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


@cache_results
def validity_check(molecule, res_molecules, res_explanations, res_confidence):
    valid_pathways = []
    valid_explanations = []
    valid_confidence = []
    for idx, smile_list in enumerate(res_molecules):
        valid = []
        if isinstance(smile_list, list):
            for smiles in smile_list:
                if is_valid_smiles(smiles):
                    if are_molecules_same(molecule, smiles):
                        logging.info(
                            f"Molecule : {molecule} is same as target molecule"
                        )
                    elif substructure_matching(smiles, molecule):
                        logging.info(
                            f"Molecule : {molecule} is substructure of target molecule"
                        )
                    else:
                        valid.append(smiles)
                else:
                    logging.info(
                        f"Molecule : {molecule} is invalid or cannot be parsed"
                    )
            if len(valid) >= 2:
                valid_pathways.append(valid)
                valid_explanations.append(res_explanations[idx])
                valid_confidence.append(res_confidence[idx])
        else:
            if is_valid_smiles(smile_list):
                if are_molecules_same(molecule, smiles):
                    logging.info("Molecule is same as target molecule")
                elif substructure_matching(smiles, molecule):
                    logging.info(
                        f"Molecule : {molecule} is substructure of target molecule {smiles}"
                    )
                else:
                    valid_pathways.append([smile_list])
                    valid_explanations.append(res_explanations[idx])
                    valid_confidence.append(res_confidence[idx])
            else:
                logging.info("Molecule is invalid or cannot be parsed")
    logging.info(
        f"Obtained {len(valid_pathways)} valid pathways after validity test: {valid_pathways}"
    )
    return valid_pathways, valid_explanations, valid_confidence


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

    # print("--------------------")
    # print(f"Storage: {storage}")
    # print(f"Fixed Dependencies: {fixed_dependencies}")
    return fixed_dependencies


def calc_yield(mol1, mol2):
    """Calculate the yield of a reaction"""
    return "#"


def calc_scalability_index(mol1, mol2):
    """Calculate the scalability index of a reaction"""
    _, type = get_reaction_type(mol1, mol2)
    return str(ENCODING_SCALABILITY[type])


def calc_chemical_formula(mol: str):
    """Calculate the chemical formula of a molecule

    Parameters
    ----------
    mol : str
        SMILES string of the molecule

    Returns
    -------
    str
        molecular formula of the molecule
    """
    try:
        formula = CalcMolFormula(Chem.MolFromSmiles(mol))
    except:
        formula = "N/A"
        logging.info(f"Error in calculating formula: {mol}")
    return formula


def calc_mol_wt(mol: str) -> float:
    """Calculate the molecular weight of a molecule

    Parameters
    ----------
    mol : str
        SMILES string of the molecule

    Returns
    -------
    float
        molecular weight of the molecule
    """
    try:
        mol_wt = ExactMolWt(Chem.MolFromSmiles(mol))
    except:
        mol_wt = 0.0
        logging.info(f"Error in calculating molecular weight: {mol}")
    return mol_wt


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


def compute_fingerprint(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,
                                                        radius,
                                                        nBits=nBits)
    return list(fingerprint)


def sub_structure_matching(target_smiles: str, query_smiles: str) -> bool:
    """Check if the query substructure is present in the target molecule"""
    target_molecule = Chem.MolFromSmiles(target_smiles)
    query_molecule = Chem.MolFromSmiles(query_smiles)

    if target_molecule.HasSubstructMatch(query_molecule):
        return True
    else:
        return False


def get_reaction_type(mol1, mol2):
    """Get the reaction type of a reaction"""
    clf = joblib.load('../reaction_prediction/rfc.pkl')
    mol1_fingerprint = compute_fingerprint(mol1)
    mol2_fingerprint = compute_fingerprint(mol2)
    reaction_type = clf.predict([mol1_fingerprint + mol2_fingerprint])
    return REACTION_ENCODING_NAMES[reaction_type[0]], reaction_type[0]


def are_molecules_same(smiles1: str, smiles2: str) -> bool:
    # Convert SMILES strings to RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        raise ValueError("Invalid SMILES string provided.")

    # Get canonical SMILES for both molecules
    canonical_smiles1 = Chem.MolToSmiles(mol1, canonical=True)
    canonical_smiles2 = Chem.MolToSmiles(mol2, canonical=True)

    # Alternatively, compare molecular fingerprints
    fingerprint1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1,
                                                                  radius=2,
                                                                  nBits=1024)
    fingerprint2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2,
                                                                  radius=2,
                                                                  nBits=1024)

    # Check if canonical SMILES or fingerprints match
    if canonical_smiles1 == canonical_smiles2:
        return True
    elif fingerprint1 == fingerprint2:
        return True
    else:
        return False


# recursive function to run pritvi
def rec_run_prithvi(molecule):
    solved, result_dict = run_az(molecule)
    result_dict = result_dict[0]
    if not solved:

        logging.info(f"AZ failed for {molecule}, running LLM")
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
        logging.info(f"LLM returned {out_pathways}")
        logging.info(f"LLM explained {out_explained}")
        for pathway in out_pathways:
            if isinstance(pathway, list):
                temp_stat = []
                for mol in pathway:
                    res, stat = rec_run_prithvi(mol)
                    if stat:
                        temp_stat.append(True)
                        result_dict['children'][0]['children'].append(res)
                logging.info(f"temp_stat: {temp_stat}")
                if all(temp_stat):
                    solved = True
            else:
                res, solved = rec_run_prithvi(pathway)
                result_dict['children'][0]['children'].append(res)
            if solved:
                logging.info('breaking')
                break
    else:
        logging.info(f"AZ solved {molecule}")
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved


def run_prithvi(molecule):
    result_dict, solved = rec_run_prithvi(molecule)
    output_data = format_output(result_dict)
    return output_data
