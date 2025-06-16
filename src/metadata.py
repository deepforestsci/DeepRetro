"""LLM-based agents for predicting chemical reaction metadata.

This module provides a suite of 'agent' functions that leverage Large Language
Models (LLMs) via the `litellm` library to predict various metadata associated
with a given chemical reaction (defined by reactants and products).

The predicted metadata categories include:
1.  **Reagents:** Likely reagents involved in the transformation.
2.  **Reaction Conditions:** Plausible temperature, pressure, solvent, and reaction time.
3.  **Literature References:** Potentially relevant scientific publications.

Each category has a pair of functions:

*   `*_agent()`: A higher-level function that takes structured input (often lists
    of dictionaries containing SMILES and other info), calls the corresponding
    `*_llm_call()` function, processes the LLM's raw output, validates it (e.g.,
    SMILES validity), and formats it into a structured dictionary or list.
    These agent functions are decorated with `@cache_results` for efficiency.
*   `*_llm_call()`: A lower-level function that constructs the specific prompt for
    the LLM (using templates from `src.variables`), makes the API call to the LLM,
    and returns the raw (but parsed from string) response.

Langfuse is integrated for logging LLM interactions.
"""

import ast
import litellm
from typing import Optional
from dotenv import load_dotenv
from litellm import completion
from src.variables import REAGENT_USER_PROMPT, REAGENT_SYS_PROMPT
from src.variables import CONDITIONS_USER_PROMPT, CONDITIONS_SYS_PROMPT
from src.variables import LITERATURE_USER_PROMPT, LITERATURE_SYS_PROMPT
from src.cache import cache_results
from src.utils.utils_molecule import calc_mol_wt, is_valid_smiles, calc_chemical_formula
from src.utils.job_context import logger as context_logger

load_dotenv()

# set the success callback to langfuse for logging
litellm.success_callback = ["langfuse"]
litellm.drop_params = True

metadata = {
    "generation_name": "prod",  # set langfuse generation name
    "project": "Retrosynthesis",  # set langfuse project name
    "version": "0.0.3",  # set langfuse version
    "trace_name": "prod",  # set langfuse Trace Name
    "trace_user_id": "sv",  # set langfuse Trace User ID
    "session_id": "metadata",  # set langfuse Session ID
}


@cache_results
def reagent_agent(reactants: list[dict],
                  product: list[dict],
                  LLM: str = "claude-3-opus-20240229",
                  temperature: float = 0.0):
    """Predicts reagents for a reaction and enriches them with metadata.

    This agent calls an LLM to suggest reagents given reactants and a product.
    It then validates the SMILES of suggested reagents and calculates their
    chemical formula and molecular weight.

    Args:
        reactants (list[dict]):
            List of reactant dictionaries. Each dictionary must have a 'smiles'
            key with the reactant's SMILES string (e.g., `[{'smiles': 'CCO'}]`).
        product (list[dict]):
            List containing a single product dictionary. The dictionary must have
            a 'smiles' key (e.g., `[{'smiles': 'CC(=O)O'}]`).
        LLM (str, optional):
            Identifier for the LLM to use. Defaults to "claude-3-opus-20240229".
        temperature (float, optional):
            Temperature for LLM sampling. Defaults to 0.0 for deterministic output.

    Returns:
        tuple[int, list[dict] | str]:
            A tuple containing a status code and the result.

            -   If successful (status 200): A list of reagent dictionaries.
                Each dictionary has the structure::

                    {
                        'smiles': '<reagent_SMILES>',
                        'reagent_metadata': {'name': '', 'chemical_formula': '...', 'mass': ...}
                    }

            -   If an error occurs (e.g., LLM call fails, parsing error, no valid
                reagents found): A status code (e.g., 404) and an empty string \"\".
    """
    logger = context_logger.get()
    product_smiles = product[0]['smiles']
    reactants_smiles = [r['smiles'] for r in reactants]
    status, res = reagent_llm_call(reactants_smiles, product_smiles, LLM,
                                   temperature)
    # logger.info(f"Received reagents from LLM: {res}")

    # Make sure the LLM call was successful
    if status != 200:
        return status, ""

    # Parse the reagents and add metadata
    try:
        reagents = res['data']
        reagent_expl = res['explanation']
    except Exception as e:
        logger.info(f"Error in parsing reagents: {e}")
        return 404, ""

    # Filter out invalid reagents
    try:
        reagents = [r for r in reagents if is_valid_smiles(r)]
    except Exception as e:
        logger.info(f"Error in filtering reagents: {e}")
        return 404, ""

    # Flag if no reagents are found
    if len(reagents) == 0:
        return 404, ""

    # Add metadata to the reagents
    res_final = []
    try:
        for _, reagent in enumerate(reagents):
            res_final.append({
                "smiles": reagent,
                "reagent_metadata": {
                    "name": "",
                    "chemical_formula": calc_chemical_formula(reagent),
                    "mass": calc_mol_wt(reagent),
                }
            })
    except Exception as e:
        logger.info(f"Error in adding metadata to reagents: {e}")
        return 404, ""
    return 200, res_final


@cache_results
def reagent_llm_call(reactants: list[str],
                     product: str,
                     LLM: str = "claude-3-opus-20240229",
                     temperature: float = 0.0):
    """Calls an LLM to predict reagents for a given reaction.

    Uses system and user prompts (`REAGENT_SYS_PROMPT`, `REAGENT_USER_PROMPT`)
    to query the LLM for suggested reagents based on reactant and product SMILES.

    Args:
        reactants (list[str]): List of SMILES strings for the reactants.
        product (str): SMILES string for the product.
        LLM (str, optional):
            Identifier for the LLM to use. Defaults to "claude-3-opus-20240229".
        temperature (float, optional):
            Temperature for LLM sampling. Defaults to 0.0.

    Returns:
        tuple[int, dict | str]:
            A tuple containing a status code and the result.

            -   If successful (status 200): A dictionary parsed from the LLM's JSON
                response. Expected format: `{'data': ['reagent1_SMILES', ...],
                'explanation': '...'}`.
            -   If an LLM call or parsing fails (status 404): An empty string "".
    """
    logger = context_logger.get()
    user_prompt = REAGENT_USER_PROMPT.replace('{product}', product)
    user_prompt = user_prompt.replace('{reactants}', ', '.join(reactants))
    messages = [{
        "role": "system",
        "content": REAGENT_SYS_PROMPT
    }, {
        "role": "user",
        "content": user_prompt
    }]
    try:
        response = completion(model=LLM,
                              messages=messages,
                              max_completion_tokens=4096,
                              temperature=temperature,
                              seed=42,
                              top_p=0.9,
                              metadata=metadata)
        res_text = response.choices[0].message.content
    except Exception as e:
        logger.info(f"Error in calling {LLM}: {e}")
        logger.info(f"Retrying call to {LLM}")
        try:
            response = completion(model=LLM,
                                  messages=messages,
                                  max_completion_tokens=4096,
                                  temperature=temperature,
                                  seed=42,
                                  top_p=0.9)
            res_text = response.choices[0].message.content
        except Exception as e:
            logger.info(f"2nd Error in calling {LLM}: {e}")
            logger.info(f"Exiting call to {LLM}")
            return 404, ""
    logger.info(f"Received response from LLM: {res_text}")
    res = ast.literal_eval(res_text)
    return 200, res


@cache_results
def conditions_agent(reactants: list[dict],
                     product: list[dict],
                     reagents: list[dict],
                     LLM: str = "claude-3-opus-20240229",
                     temperature: float = 0.0):
    """Predicts reaction conditions using an LLM.

    Given reactants, product, and reagents, this agent calls an LLM to suggest
    reaction conditions (temperature, pressure, solvent, time).

    Args:
        reactants (list[dict]):
            List of reactant dictionaries, each with a 'smiles' key.
        product (list[dict]):
            List with a single product dictionary, with a 'smiles' key.
        reagents (list[dict]):
            List of reagent dictionaries, each with a 'smiles' key.
        LLM (str, optional):
            Identifier for the LLM. Defaults to "claude-3-opus-20240229".
        temperature (float, optional):
            LLM sampling temperature. Defaults to 0.0.

    Returns:
        tuple[int, dict | str]:
            A tuple containing a status code and the result.

            -   If successful (status 200): A dictionary containing the predicted
                conditions, e.g., `{'temperature': '25 C', 'pressure': '1 atm',
                'solvent': 'Water', 'time': '2 h'}`.
            -   If an error occurs (LLM call, parsing): Status code (e.g., 404)
                and an empty string "".
    """
    logger = context_logger.get()
    product_smiles = product[0]['smiles']
    reactants_smiles = [r['smiles'] for r in reactants]
    reagents_smiles = [r['smiles'] for r in reagents]
    status, res = conditions_llm_call(reactants_smiles, product_smiles,
                                      reagents_smiles, LLM, temperature)
    # logger.info(f"Received conditions from LLM: {res}")

    # Make sure the LLM call was successful
    if status != 200:
        return status, ""

    # Make sure the conditions are present
    try:
        temp = res['temperature']
        pressure = res['pressure']
        solvent = res['solvent']
        time = res['time']
    except Exception as e:
        logger.info(f"Error in parsing conditions: {e}")
        return 404, ""

    return 200, res


@cache_results
def conditions_llm_call(reactants: list[str],
                        product: str,
                        reagents: list[str],
                        LLM: str = "claude-3-opus-20240229",
                        temperature: float = 0.0):
    """Calls an LLM to predict reaction conditions.

    Uses `CONDITIONS_SYS_PROMPT` and `CONDITIONS_USER_PROMPT` to query the LLM.

    Args:
        reactants (list[str]): List of SMILES strings for reactants.
        product (str): SMILES string for the product.
        reagents (list[str]): List of SMILES strings for reagents.
        LLM (str, optional):
            Identifier for the LLM. Defaults to "claude-3-opus-20240229".
        temperature (float, optional):
            LLM sampling temperature. Defaults to 0.0.

    Returns:
        tuple[int, dict | str]:
            A tuple containing a status code and the result.

            -   If successful (status 200): A dictionary parsed from the LLM's JSON
                response, e.g., `{'temperature': '25 C', ...}`.
            -   If an LLM call or parsing fails (status 404): An empty string "".
    """
    logger = context_logger.get()
    user_prompt = CONDITIONS_USER_PROMPT.replace('{product}', product)
    user_prompt = user_prompt.replace('{reactants}', ', '.join(reactants))
    user_prompt = user_prompt.replace('{reagents}', ', '.join(reagents))
    messages = [{
        "role": "system",
        "content": CONDITIONS_SYS_PROMPT
    }, {
        "role": "user",
        "content": user_prompt
    }]
    try:
        response = completion(model=LLM,
                              messages=messages,
                              max_completion_tokens=4096,
                              temperature=temperature,
                              seed=42,
                              top_p=0.9,
                              metadata=metadata)
        res_text = response.choices[0].message.content
    except Exception as e:
        logger.info(f"Error in calling {LLM}: {e}")
        logger.info(f"Retrying call to {LLM}")
        try:
            response = completion(model=LLM,
                                  messages=messages,
                                  max_completion_tokens=4096,
                                  temperature=temperature,
                                  seed=42,
                                  top_p=0.9)
            res_text = response.choices[0].message.content
        except Exception as e:
            logger.info(f"2nd Error in calling {LLM}: {e}")
            logger.info(f"Exiting call to {LLM}")
            return 404, ""
    logger.info(f"Received response from LLM: {res_text}")
    res_text = ast.literal_eval(res_text)
    return 200, res_text


@cache_results
def literature_agent(reactants: list[dict],
                     product: list[dict],
                     reagents: list[dict],
                     conditions: dict,
                     LLM: str = "claude-3-opus-20240229",
                     temperature: float = 0.0):
    """Predicts literature references for a reaction using an LLM.

    This agent calls an LLM to suggest relevant literature references given
    reactants, product, reagents, and reaction conditions.
    It directly handles the LLM interaction using `LITERATURE_SYS_PROMPT`
    and `LITERATURE_USER_PROMPT`.

    Args:
        reactants (list[dict]):
            List of reactant dictionaries, each with a 'smiles' key.
        product (list[dict]):
            List with a single product dictionary, with a 'smiles' key.
        reagents (list[dict]):
            List of reagent dictionaries, each with a 'smiles' key.
        conditions (dict):
            A dictionary describing reaction conditions, typically the output
            from `conditions_agent` (e.g., `{'temperature': '25 C', ...}`).
            This is converted to a string for the prompt.
        LLM (str, optional):
            Identifier for the LLM. Defaults to "claude-3-opus-20240229".
        temperature (float, optional):
            LLM sampling temperature. Defaults to 0.0.

    Returns:
        tuple[int, list | dict | str]:
            A tuple containing a status code and the result.

            -   If successful (status 200): A list or dictionary containing the
                predicted literature references. The exact structure depends on the
                LLM prompt, but typically includes fields like DOI, title, authors,
                journal, year for each reference (e.g., `[{'doi': '...', ...}]`).
            -   If an error occurs (LLM call, parsing): Status code (e.g., 404)
                and an empty string "".
    """
    logger = context_logger.get()
    product_smiles = product[0]['smiles']
    reactants_smiles = [r['smiles'] for r in reactants]
    reagents_smiles = [r['smiles'] for r in reagents]
    user_prompt = LITERATURE_USER_PROMPT.replace('{product}', product_smiles)
    user_prompt = user_prompt.replace('{reactants}',
                                      ', '.join(reactants_smiles))
    user_prompt = user_prompt.replace('{reagents}', ', '.join(reagents_smiles))
    user_prompt = user_prompt.replace('{conditions}', str(conditions))
    messages = [{
        "role": "system",
        "content": LITERATURE_SYS_PROMPT
    }, {
        "role": "user",
        "content": user_prompt
    }]
    try:
        response = completion(model=LLM,
                              messages=messages,
                              max_completion_tokens=4096,
                              temperature=temperature,
                              seed=42,
                              top_p=0.9,
                              metadata=metadata)
        res_text = response.choices[0].message.content
    except Exception as e:
        logger.info(f"Error in calling {LLM}: {e}")
        logger.info(f"Retrying call to {LLM}")
        try:
            response = completion(model=LLM,
                                  messages=messages,
                                  max_completion_tokens=4096,
                                  temperature=temperature,
                                  seed=42,
                                  top_p=0.9)
            res_text = response.choices[0].message.content
        except Exception as e:
            logger.info(f"2nd Error in calling {LLM}: {e}")
            logger.info(f"Exiting call to {LLM}")
            return 404, ""
    logger.info(f"Received response from LLM: {res_text}")

    # convert the response to a dictionary
    try:
        res_text = ast.literal_eval(res_text)
    except Exception as e:
        logger.info(f"Error in parsing literature output: {e}")
        return 404, ""

    # Parse the literature reaction
    try:
        res_lit = res_text['literature_reaction']
    except Exception as e:
        logger.info(f"Error in parsing literature reaction: {e}")
        return 404, ""
    return 200, res_lit
