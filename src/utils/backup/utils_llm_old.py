"""Utils for calling LLM model and processing the results"""
import ast
import os
import logging
from typing import Optional
from anthropic import Anthropic
from dotenv import load_dotenv
from src.variables import USER_PROMPT_OLD, SYS_PROMPT_OLD
from src.cache import cache_results
from src.utils.utils_molecule import calc_mol_wt, validity_check

# load environment variables
load_dotenv()
ANTHROPIC_API_KEY = os.getenv('ANTHROPIC_API_KEY')

logger = logging.getLogger(__name__)

client = Anthropic(api_key=ANTHROPIC_API_KEY)


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
                USER_PROMPT_OLD.replace('{target_smiles}', molecule)
            }]
        }]
    try:
        message = client.messages.create(model=LLM,
                                         max_tokens=1024,
                                         temperature=temperature,
                                         system=SYS_PROMPT_OLD,
                                         messages=messages,
                                         top_p=0.9)
        res_text = message.content[0].text
    except Exception as e:
        logger.info(f"Error in calling LLM: {e}")
        message = client.messages.create(model=LLM,
                                         max_tokens=1024,
                                         temperature=temperature,
                                         system=SYS_PROMPT_OLD,
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

    output_pathways.sort(key=lambda x: calc_mol_wt(x[0]))
    # TODO: fix the matching of output pathways and confidence
    return output_pathways, output_explanations, output_confidence
