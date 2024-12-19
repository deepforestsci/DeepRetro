"""Code for Metadata agents
Metadata covers the following:
1. Nearest Literature (1-2) - 1-2 nearest literature references, to contain the doi, title, authors, journal, year
2. Reagents (1-2) - 1-2 reagents used in the reaction, to contain the SMILES
3. Reaction conditions - Conditions used in the reaction, to contain the temperature, pressure, solvent, time

Each of these will be a separate agent, and will be called by the main agent to get the metadata for the reaction
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
from src.utils.utils_molecule import calc_mol_wt, validity_check
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
def reagent_agent(reactants: list[str],
                  product: str,
                  LLM: str = "claude-3-opus-20240229",
                  temperature: float = 0.0):
    """Calls the LLM model to predict the reagents used in the reaction

    Parameters
    ----------
    reactants : list[str]
        List of reactants SMILES
    product : str
        Product SMILES
    LLM : str, optional
        LLM model to use, by default "claude-3-opus-20240229"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, list[str]
        Status code and list of reagents SMILES
    """
    logger = context_logger.get()
    user_prompt = REAGENT_USER_PROMPT.replace('{target_smiles}', product)
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
def conditions_agent(reactants: list[str],
                     product: str,
                     reagents: list[str],
                     LLM: str = "claude-3-opus-20240229",
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
        LLM model to use, by default "claude-3-opus-20240229"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, str
        Status code and response text
    """
    logger = context_logger.get()
    user_prompt = CONDITIONS_USER_PROMPT.replace('{target_smiles}', product)
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
    return 200, res_text


@cache_results
def literature_agent(reactants: list[str],
                     product: str,
                     reagents: list[str],
                     conditions: str,
                     LLM: str = "claude-3-opus-20240229",
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
        LLM model to use, by default "claude-3-opus-20240229"
    temperature : float, optional
        Temperature for the LLM model, by default 0.0

    Returns
    -------
    int, str
        Status code and response text
    """
    logger = context_logger.get()
    user_prompt = LITERATURE_USER_PROMPT.replace('{target_smiles}', product)
    user_prompt = user_prompt.replace('{reactants}', ', '.join(reactants))
    user_prompt = user_prompt.replace('{reagents}', ', '.join(reagents))
    user_prompt = user_prompt.replace('{conditions}', conditions)
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
    return 200, res_text
