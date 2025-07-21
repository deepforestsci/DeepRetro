"""
Code for Metadata agents

Metadata covers:
1. Nearest Literature (1-2) - 1-2 nearest literature references, to contain the doi, title, authors, journal, year
2. Reagents (1-2) - 1-2 reagents used in the reaction, to contain the SMILES
3. Reaction conditions - Conditions used in the reaction, to contain the temperature, pressure, solvent, time

Each of these will be a separate agent, and will be called by the main agent to get the metadata for the reaction.
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
                  LLM: str = "claude-opus-4-20250514",
                  temperature: float = 0.0):
    """
    Calls the LLM model to predict the reagents used in the reaction.

    Parameters
    ----------
    reactants : list[dict]
        List of reactants dict with SMILES and metadata
    product : list[dict]
        Product dict with SMILES and metadata
    LLM : str, optional
        LLM model to use, by default "claude-opus-4-20250514"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, list[str]
        Status code and list of reagents SMILES

    Examples
    --------
    >>> reactants = [{'smiles': 'CCO'}]
    >>> product = [{'smiles': 'CC=O'}]
    >>> status, reagents = reagent_agent(reactants, product)  # doctest: +SKIP
    >>> status  # doctest: +SKIP
    200
    >>> isinstance(reagents, list)  # doctest: +SKIP
    True
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
                     LLM: str = "claude-opus-4-20250514",
                     temperature: float = 0.0):
    """Calls the LLM model to predict the reagents used in the reaction

    Parameters
    ----------
    reactants : list[str]
        List of reactants SMILES
    product : list[str]
        Product SMILES
    LLM : str, optional
        LLM model to use, by default "claude-opus-4-20250514"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, list[str]
        Status code and list of reagents SMILES
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
                     LLM: str = "claude-opus-4-20250514",
                     temperature: float = 0.0):
    """Calls the LLM model to predict the reaction conditions

    Parameters
    ----------
    reactants : list[str]
        List of reactants SMILES
    product : str
        Product SMILES
    reagents : list[str]
        List of reagents SMILES
    LLM : str, optional
        LLM model to use, by default "claude-opus-4-20250514"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, str
        Status code and response text
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
                        LLM: str = "claude-opus-4-20250514",
                        temperature: float = 0.0):
    """Calls the LLM model to predict the reaction conditions

    Parameters
    ----------
    reactants : list[str]
        List of reactants SMILES
    product : str
        Product SMILES
    reagents : list[str]
        List of reagents SMILES
    LLM : str, optional
        LLM model to use, by default "claude-opus-4-20250514"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, str
        Status code and response text
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
def literature_agent(reactants: list[str],
                     product: str,
                     reagents: list[str],
                     conditions: str,
                     LLM: str = "claude-opus-4-20250514",
                     temperature: float = 0.0):
    """Calls the LLM model to predict the nearest literature references

    Parameters
    ----------
    reactants : list[str]
        List of reactants SMILES
    product : str
        Product SMILES
    reagents : list[str]
        List of reagents SMILES
    conditions : str
        Reaction conditions
    LLM : str, optional
        LLM model to use, by default "claude-opus-4-20250514"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, str
        Status code and response text
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
