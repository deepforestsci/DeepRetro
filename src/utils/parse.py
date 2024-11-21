""" This module contains functions to parse the input data to extract steps and their dependencies. """
from typing import Dict, List, Optional
from src.utils.utils_molecule import validity_check, calc_chemical_formula, calc_mol_wt
from src.utils.utils_molecule import calc_confidence_estimate, calc_scalability_index
from src.variables import BASIC_MOLECULES, ENCODING_SCALABILITY


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
                dependency_list[str(parent_id)].append(str(step_id))
            else:
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

        step_list[parent_id - 1]["reactionmetrics"][0][
            "scalabilityindex"] = calc_scalability_index(
                data['smiles'],
                step_list[parent_id - 1]["products"][0]["smiles"])

    if step is not None:
        step_list.append(step)
        if 'children' in data:
            for child in data['children']:
                if 'children' in child:
                    for c in child['children']:
                        parse_step(c, include_metadata, step_list,
                                   dependency_list, step_id)
                else:
                    parse_step(child, include_metadata, step_list,
                               dependency_list, step_id)

    return {"dependencies": dependency_list, "steps": step_list}


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


def format_output(data):
    """Format the output data for visualization"""
    output_data = parse_step(data)
    # fix_dependencies(output_data['dependencies'], output_data['steps'])
    output_data['dependencies'] = fix_dependencies(output_data['dependencies'],
                                                   output_data['steps'])
    return output_data
