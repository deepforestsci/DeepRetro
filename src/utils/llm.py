import ast
import litellm
from typing import Optional
from dotenv import load_dotenv
from litellm import completion
from src.variables import USER_PROMPT, SYS_PROMPT
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
    "session_id": "prod",  # set langfuse Session ID
}


@cache_results
def call_LLM(molecule: str,
             LLM: str = "claude-3-opus-20240229",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None):
    """Calls the LLM model to predict the next step"""
    logger = context_logger.get()
    # logger.info(f"Calling {LLM} with molecule: {molecule}")
    if messages is None:
        messages = [{
            "role": "system",
            "content": SYS_PROMPT
        }, {
            "role": "user",
            "content": USER_PROMPT.replace('{target_smiles}', molecule)
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


def split_cot_json(res_text: str) -> tuple[int, list[str], str]:
    """Parse the LLM response to extract the thinking steps and json content

    Parameters
    ----------
    res_text : str
        The response text from the LLM model

    Returns
    -------
    tuple[int, list[str], str]
        The status code, thinking steps and json content
    """
    logger = context_logger.get()
    try:
        # extract the content within <cot> </cot> tags as thinking content
        thinking_content = res_text[res_text.find("<cot>\n") +
                                    6:res_text.find("</cot>")]
        # split the thinking content into individual steps based on the <thinking> </thinking> tags
        thinking_steps = thinking_content.split("<thinking>\n")[1:]
        thinking_steps = [
            step[:step.find("</thinking>")] for step in thinking_steps
        ]
        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
    except Exception as e:
        logger.info(f"Error in parsing LLM response: {e}")
        return 500, [], ""
    return 200, thinking_steps, json_content


def validate_split_json(
        json_content: str) -> tuple[int, list[str], list[str], list[int]]:
    """Validate the split json content from LLM response

    Parameters
    ----------
    json_content : str
        The json content from the LLM response

    Returns
    -------
    tuple[int, list[str], list[str], list[int]]
        The status code, list of molecules, list of explanations and list of confidence scores
    """
    logger = context_logger.get()
    try:
        result_list = ast.literal_eval(json_content)
        res_molecules = result_list['data']
        res_explanations = result_list['explanation']
        res_confidence = result_list['confidence_scores']
    except Exception as e:
        logger.info(f"Error in parsing response: {e}")
        return 500, [], [], []
    return 200, res_molecules, res_explanations, res_confidence


def llm_pipeline(
    molecule: str,
    LLM: str = "claude-3-opus-20240229",
    messages: Optional[list[dict]] = None
) -> tuple[list[str], list[str], list[float]]:
    """Pipeline to call LLM and validate the results

    Parameters
    ----------
    molecule : str
        The target molecule for retrosynthesis
    LLM : str, optional
        LLM to be used for retrosynthesis , by default "claude-3-opus-20240229"
    messages : Optional[list[dict]], optional
        Conversation history, by default None

    Returns
    -------
    tuple[list[str], list[str], list[float]]
        The output pathways, explanations and confidence scores
    """
    logger = context_logger.get()
    output_pathways = []
    run = 0.0
    while (output_pathways == [] and run < 0.6):
        logger.info(f"Calling LLM with molecule: {molecule} and run: {run}")
        status_code, res_text = call_LLM(molecule,
                                         LLM,
                                         messages=messages,
                                         temperature=run)
        if status_code == 200:
            status_code, thinking_steps, json_content = split_cot_json(
                res_text)
            if status_code == 200:
                status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
                    json_content)
                if status_code == 200:
                    output_pathways, output_explanations, output_confidence = validity_check(
                        molecule, res_molecules, res_explanations,
                        res_confidence)
                    logger.info(f"Output Pathways: {output_pathways},\n\
                            Output Explanations: {output_explanations},\n\
                                Output Confidence: {output_confidence}")
                    run += 0.1
                else:
                    logger.info(
                        f"Error in validating split json content: {res_text}")
                    continue
            else:
                logger.info(f"Error in splitting cot json: {res_text}")
                continue
        else:
            logger.info(f"Error in calling LLM: {res_text}")
            continue
    return output_pathways, output_explanations, output_confidence
