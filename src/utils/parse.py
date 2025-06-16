"""Parses hierarchical retrosynthetic route data into a flat list of steps and dependencies.

This module is designed to process data representing retrosynthetic pathways,
which are typically structured as a tree or graph where nodes are molecules
and reactions. The primary function, `parse_step`, recursively traverses this
input data (expected to have 'smiles' and 'children' keys) to extract individual
reaction steps and identify their interdependencies.

The output is a dictionary containing a list of 'steps' and a 'dependencies'
map. Each step details reactants, products, and associated metadata (like
chemical formula, molecular weight, and preliminary reaction metrics).
The `fix_dependencies` function further refines the dependency information to a
format suitable for downstream processing or visualization, mapping each step to
the steps that produce its reactants.

This module relies on utility functions from `src.utils.utils_molecule` for
calculating chemical properties and `src.variables` for predefined basic
molecules and scalability encodings.
"""
from typing import Dict, List, Optional
from src.utils.utils_molecule import validity_check, calc_chemical_formula, calc_mol_wt
from src.utils.utils_molecule import calc_confidence_estimate, calc_scalability_index
from src.variables import BASIC_MOLECULES, ENCODING_SCALABILITY


def parse_step(data,
               include_metadata: Optional[bool] = None,
               step_list: Optional[List] = None,
               dependency_list: Optional[Dict] = None,
               parent_id: Optional[int] = None) -> Dict:
    """Recursively parses hierarchical route data to extract reaction steps and dependencies.

    This function traverses a nested dictionary structure (typically from AiZynthFinder
    or a similar source) representing a retrosynthetic tree. It flattens this
    tree into a list of reaction steps and a dictionary defining the dependencies
    between these steps (parent step producing a product which is a reactant in a child step).

    Args:
        data (Dict):
            A dictionary representing a node in the retrosynthetic tree. Expected keys
            include 'smiles' (str, SMILES of the molecule) and optionally 'children'
            (List[Dict], child nodes which can be reactants or intermediate reactions).
            Other keys like 'is_reaction' and 'metadata' (containing, e.g.,
            'policy_probability') are also processed if present.
        include_metadata (Optional[bool], optional):
            A flag to determine if metadata (chemical formula, molecular weight) should be
            calculated and included for products, reactants, and reagents. Defaults to True.
        step_list (Optional[List], optional):
            Accumulator list for the extracted reaction steps. This is used internally
            during recursion. Users typically call this function with `step_list=None`.
            Defaults to None, initializing an empty list.
        dependency_list (Optional[Dict], optional):
            Accumulator dictionary for parent-child step dependencies. This is used
            internally. Users typically call with `dependency_list=None`.
            Defaults to None, initializing an empty dictionary.
        parent_id (Optional[int], optional):
            The ID of the parent step in the recursion. Used internally to build the
            dependency list. Defaults to None for the root call.

    Returns:
        Dict:
            A dictionary containing two keys:
            - 'dependencies' (Dict[str, List[str]]): Maps a parent step ID (str) to a
              list of its child step IDs (str). This represents the initial hierarchical
              dependency.
            - 'steps' (List[Dict]): A list where each dictionary represents a reaction step.
              The structure of each step dictionary is:
                - 'step' (str): A unique identifier for the step (e.g., "1", "2").
                - 'reactants' (List[Dict]): List of reactants. Each reactant dict contains:
                    - 'smiles' (str): SMILES string.
                    - 'reactant_metadata' (Dict): Contains 'name' (str, currently empty),
                      'chemical_formula' (str), 'mass' (float).
                - 'reagents' (List[Dict]): List of reagents (molecules from `BASIC_MOLECULES`).
                  Structure similar to reactants, with 'reagent_metadata'.
                - 'products' (List[Dict]): List of products (usually one per step from input).
                  Structure similar to reactants, with 'product_metadata'.
                - 'conditions' (List): Currently an empty list. Intended for reaction conditions.
                - 'reactionmetrics' (List[Dict]): List containing one dictionary with reaction
                  metrics:
                    - 'scalabilityindex' (str): Calculated scalability index.
                    - 'confidenceestimate' (float): Calculated from policy probability.
                    - 'closestliterature' (str): Currently empty.
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
    """Transforms the dependency list to map steps to their precursor steps.

    The initial `dependencies` map (from `parse_step`) shows a parent step ID
    pointing to its children step IDs. This function inverts this relationship
    for a more direct graph representation: for each step, it identifies which
    other steps produce its reactants.

    Args:
        dependencies (Dict[str, List[str]]):
            The initial dependency map. While passed as an argument, this specific
            implementation does not directly use this input `dependencies` dictionary.
            It re-calculates dependencies based on the `step_list`.
            (Note: The original `dependencies` argument might be from a previous
            design or intended for a different type of fixing. The current code
            derives dependencies solely from `step_list`.)
        step_list (List[Dict]):
            A list of step dictionaries, as produced by `parse_step`. Each step
            dictionary must contain 'products' with 'smiles' and a 'step' ID, and
            'reactants' with 'smiles'.

    Returns:
        Dict[str, List[str]]:
            A dictionary where each key is a step ID (str), and the value is a list
            of step IDs (str) that produce the reactants for the key step.
            Effectively, `fixed_dependencies[step_X_id] = [step_A_id, step_B_id]` means
            step A and step B produce reactants needed for step X.
    """
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
    """Formats raw hierarchical route data into a structured step and dependency format.

    This function serves as the main entry point for parsing. It first calls
    `parse_step` to convert the input hierarchical data into a list of steps
    and an initial parent-to-child dependency map. Then, it calls
    `fix_dependencies` to transform this dependency map into a format where
    each step is mapped to the steps that produce its reactants.

    Args:
        data (Dict):
            The raw input data representing a retrosynthetic route, typically a
            nested dictionary structure with 'smiles' and 'children' keys at each
            node (e.g., output from AiZynthFinder or a similar tool).

    Returns:
        Dict:
            A dictionary with two keys:
            - 'dependencies' (Dict[str, List[str]]): The transformed dependency map
              from `fix_dependencies`, where each step ID maps to a list of IDs of
              steps that produce its reactants.
            - 'steps' (List[Dict]): The list of parsed step dictionaries from `parse_step`.
    """
    output_data = parse_step(data)
    # fix_dependencies(output_data['dependencies'], output_data['steps'])
    output_data['dependencies'] = fix_dependencies(output_data['dependencies'],
                                                   output_data['steps'])
    return output_data
