import os
import ast
import litellm
from typing import Optional
from dotenv import load_dotenv
from litellm import completion
from src.variables import OPENAI_MODELS, DEEPSEEK_MODELS
from src.variables import USER_PROMPT_V3 as USER_PROMPT, SYS_PROMPT_V3 as SYS_PROMPT
from src.variables import USER_PROMPT_OPENAI, SYS_PROMPT_OPENAI
from src.variables import USER_PROMPT_DEEPSEEK, SYS_PROMPT_DEEPSEEK
from src.variables import ADDON_PROMPT_7_MEMBER
from src.cache import cache_results
from src.utils.utils_molecule import validity_check, detect_seven_member_rings
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
ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING",
                                    "true").lower() == "false" else True


def log_message(message: str, logger=None):
    """Log the message"""
    if logger is not None:
        logger.info(message)
    else:
        print(message)


@cache_results
def call_LLM(molecule: str,
             LLM: str = "claude-3-opus-20240229",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None):
    """Calls the LLM model to predict the next step"""
    logger = context_logger.get() if ENABLE_LOGGING else None
    log_message(f"Calling {LLM} with molecule: {molecule}", logger)

    if detect_seven_member_rings(molecule):
        log_message(f"Detected seven member ring in molecule: {molecule}",
                    logger)
        add_on = ADDON_PROMPT_7_MEMBER
    else:
        add_on = ""
    add_on = ""

    if LLM in DEEPSEEK_MODELS:
        if messages is None:
            messages = [{
                "role": "system",
                "content": SYS_PROMPT_DEEPSEEK + add_on
            }, {
                "role":
                "user",
                "content":
                USER_PROMPT_DEEPSEEK.replace('{target_smiles}', molecule)
            }]
        max_completion_tokens = 8192
    elif LLM in OPENAI_MODELS:
        if messages is None:
            messages = [{
                "role":
                "user",
                "content":
                USER_PROMPT_OPENAI.replace('{target_smiles}', molecule)
            }]
        max_completion_tokens = 8192
    else:
        if messages is None:
            messages = [{
                "role": "system",
                "content": SYS_PROMPT + add_on
            }, {
                "role":
                "user",
                "content":
                USER_PROMPT.replace('{target_smiles}', molecule)
            }]
        max_completion_tokens = 4096
    try:
        response = completion(model=LLM,
                              messages=messages,
                              max_completion_tokens=max_completion_tokens,
                              temperature=temperature,
                              seed=42,
                              top_p=0.9,
                              metadata=metadata)
        res_text = response.choices[0].message.content
    except Exception as e:
        log_message(f"Error in calling {LLM}: {e}", logger)
        log_message(f"Retrying call to {LLM}", logger)
        try:
            response = completion(model=LLM,
                                  messages=messages,
                                  max_completion_tokens=4096,
                                  temperature=temperature,
                                  seed=42,
                                  top_p=0.9)
            res_text = response.choices[0].message.content
        except Exception as e:
            log_message(f"2nd Error in calling {LLM}: {e}", logger)
            log_message(f"Exiting call to {LLM}", logger)
            return 404, ""
    log_message(f"Received response from LLM: {res_text}", logger)
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
    logger = context_logger.get() if ENABLE_LOGGING else None
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
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 500, [], ""
    return 200, thinking_steps, json_content


def split_json_openAI(res_text: str) -> tuple[int, list[str]]:
    """Split the response text from OpenAI models to extract the molecules
    Note: OpenAI O-series models do not provide Chain of Thoughts (COT) in the response

    Parameters
    ----------
    res_text : str
        The response text from the OpenAI model

    Returns
    -------
    tuple[int, str]
        the status code and json content
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 500, [], ""
    return 200, json_content


def split_json_deepseek(res_text: str) -> tuple[int, list[str], str]:
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
    logger = context_logger.get() if ENABLE_LOGGING else None

    try:
        # extract the content within <cot> </cot> tags as thinking content
        thinking_content = res_text[res_text.find("<think>\n") +
                                    6:res_text.find("</think>")]

        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 500, [], ""
    return 200, thinking_content, json_content


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
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        result_list = ast.literal_eval(json_content)
        res_molecules = result_list['data']
        res_explanations = result_list['explanation']
        res_confidence = result_list['confidence_scores']
    except Exception as e:
        log_message(f"Error in parsing response: {e}", logger)
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
    logger = context_logger.get() if ENABLE_LOGGING else None
    output_pathways = []
    run = 0.0
    while (output_pathways == [] and run < 0.6):
        log_message(f"Calling LLM with molecule: {molecule} and run: {run}",
                    logger)
        if LLM in DEEPSEEK_MODELS and run == 0.0:
            status_code, res_text = call_LLM(molecule,
                                             LLM,
                                             messages=messages,
                                             temperature=0.5)
        elif LLM in DEEPSEEK_MODELS:
            status_code, res_text = call_LLM(molecule,
                                             "claude-3-opus-20240229",
                                             messages=messages,
                                             temperature=run)
        else:
            status_code, res_text = call_LLM(molecule,
                                             LLM,
                                             messages=messages,
                                             temperature=run)
        if status_code == 200:
            if LLM in OPENAI_MODELS:
                status_code, json_content = split_json_openAI(res_text)
            elif LLM in DEEPSEEK_MODELS:
                status_code, thinking_steps, json_content = split_json_deepseek(
                    res_text)
            else:
                status_code, thinking_steps, json_content = split_cot_json(
                    res_text)
            if status_code == 200:
                status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
                    json_content)
                if status_code == 200:
                    output_pathways, output_explanations, output_confidence = validity_check(
                        molecule, res_molecules, res_explanations,
                        res_confidence)
                    log_message(
                        f"Output Pathways: {output_pathways},\n\
                            Output Explanations: {output_explanations},\n\
                                Output Confidence: {output_confidence}",
                        logger)
                    run += 0.1
                else:
                    log_message(
                        f"Error in validating split json content: {res_text}",
                        logger)
                    continue
            else:
                log_message(f"Error in splitting cot json: {res_text}", logger)
                print(f"Error in splitting cot json: {res_text}")
                continue
        else:
            log_message(f"Error in calling LLM: {res_text}")
            continue
    return output_pathways, output_explanations, output_confidence
