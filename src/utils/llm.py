"""Utilities for Interacting with Large Language Models (LLMs).

This module centralizes all functionality related to Large Language Model (LLM)
operations within the DeepRetro application. It leverages the `litellm` library
to provide a consistent interface across various LLM providers and models.

Key responsibilities include:
- **Prompt Management (`obtain_prompt`):
  Dynamically selecting appropriate system and user prompts based on the specified
  LLM model identifier. Supports standard and "advanced" prompt versions (indicated
  by an ":adv" suffix on the LLM model string).

- **LLM API Calls (`call_LLM`):
  Encapsulating the logic for making requests to LLMs. This includes constructing
  the message payload, handling model-specific parameters (e.g., for DeepSeek or
  specific Claude versions), incorporating conditional prompt additions (like for
  7-membered rings), and implementing retry mechanisms. LLM calls are cached.

- **Response Parsing (various `split_*` functions):
  Extracting structured data from potentially complex LLM responses. This includes:
    - Separating Chain-of-Thought (COT) narratives from JSON payloads
      (e.g., `split_cot_json`).
    - Handling model-specific output formats (e.g., `split_json_openAI`,
      `split_json_deepseek`).
    - A master parsing function (`split_json_master`) to dispatch to the correct
      parser based on the model.

- **Output Validation (`validate_split_json`):
  Ensuring the integrity and correctness of the data extracted from LLM responses,
  such as validating SMILES strings and expected JSON structures.

- **Main Pipeline (`llm_pipeline`):
  Orchestrating the sequence of calling an LLM, parsing its response, validating
  the output, and potentially performing other checks like stability or hallucination
  analysis.

- **Error Handling (`get_error_log`):
  Mapping status codes to descriptive error messages.

Langfuse is integrated for logging LLM interactions, controlled by an environment
variable `ENABLE_LOGGING`.
"""
import os
import ast
import litellm
from typing import Optional
import structlog
from dotenv import load_dotenv
from litellm import completion
from src.variables import OPENAI_MODELS, DEEPSEEK_MODELS
from src.variables import USER_PROMPT, SYS_PROMPT
from src.variables import USER_PROMPT_V4, SYS_PROMPT_V4
from src.variables import USER_PROMPT_OPENAI, SYS_PROMPT_OPENAI
from src.variables import USER_PROMPT_DEEPSEEK, SYS_PROMPT_DEEPSEEK
from src.variables import ADDON_PROMPT_7_MEMBER, USER_PROMPT_DEEPSEEK_V4
from src.variables import ERROR_MAP
from src.cache import cache_results
from src.utils.utils_molecule import validity_check, detect_seven_member_rings
from src.utils.job_context import logger as context_logger
from src.utils.stability_checks import stability_checker
from src.utils.hallucination_checks import hallucination_checker

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


def log_message(message: str, logger: Optional[structlog.stdlib.BoundLogger] = None):
    """Logs a message using a provided logger or defaults to print.

    Args:
        message (str): The message to be logged.
        logger (Optional[structlog.stdlib.BoundLogger], optional):
            A structlog bound logger instance. If None, `print()` is used.
            Defaults to None.
    """
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def obtain_prompt(LLM: str) -> tuple[str, str, int]:
    """Determines the system prompt, user prompt, and max completion tokens based on the LLM model.

    The `LLM` model string can include an optional ":adv" suffix (e.g.,
    "claude-3-opus-20240229:adv") to select an "advanced" version of prompts.
    Different base prompts and token limits are defined for OpenAI, DeepSeek,
    and other (default Claude) models.

    Args:
        LLM (str):
            The LLM model identifier, optionally with an ":adv" suffix.

    Returns:
        tuple[str, str, int]:
            - `sys_prompt_final` (str): The selected system prompt.
            - `user_prompt_final` (str): The selected user prompt template.
            - `max_completion_tokens` (int): The maximum number of tokens for
              the LLM completion.
    """
    advanced_prompt = False
    detector = LLM.split(":")
    if len(detector) > 1 and detector[1] == "adv":
        advanced_prompt = True
    print(f"Advanced Prompt: {advanced_prompt}")
    if advanced_prompt:
        if LLM in DEEPSEEK_MODELS:
            sys_prompt_final = SYS_PROMPT_V4
            user_prompt_final = USER_PROMPT_DEEPSEEK_V4
            max_completion_tokens = 8192 * 2
        elif LLM in OPENAI_MODELS:
            sys_prompt_final = SYS_PROMPT_OPENAI
            user_prompt_final = USER_PROMPT_OPENAI
            max_completion_tokens = 8192
        else:
            sys_prompt_final = SYS_PROMPT_V4
            user_prompt_final = USER_PROMPT_V4
            max_completion_tokens = 4096
    else:
        if LLM in DEEPSEEK_MODELS:
            sys_prompt_final = SYS_PROMPT_DEEPSEEK
            user_prompt_final = USER_PROMPT_DEEPSEEK
            max_completion_tokens = 8192 * 2
        elif LLM in OPENAI_MODELS:
            sys_prompt_final = SYS_PROMPT_OPENAI
            user_prompt_final = USER_PROMPT_OPENAI
            max_completion_tokens = 8192
        else:
            sys_prompt_final = SYS_PROMPT
            user_prompt_final = USER_PROMPT
            max_completion_tokens = 4096
    return sys_prompt_final, user_prompt_final, max_completion_tokens


@cache_results
def call_LLM(molecule: str,
             LLM: str = "claude-3-opus-20240229",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None) -> tuple[int, str]:
    """Calls the specified LLM to predict the next retrosynthetic step(s).

    Constructs prompts using `obtain_prompt`, adds specific instructions if a
    7-membered ring is detected (for some models), and adjusts API parameters for
    certain LLM types (e.g., DeepSeek, specific Claude versions like "3-7").
    Includes a retry mechanism for the API call.
    The call is cached via `@cache_results`.

    Args:
        molecule (str):
            The SMILES string of the target molecule for retrosynthesis.
        LLM (str, optional):
            The LLM model identifier (e.g., "claude-3-opus-20240229", or
            "model-name:adv" for advanced prompts). Defaults to
            "claude-3-opus-20240229". The base model name is extracted before the call.
        temperature (float, optional):
            The temperature for LLM sampling. Defaults to 0.0.
        messages (Optional[list[dict]], optional):
            A list of message dictionaries for the conversation history, each like
            `{"role": "user/system/assistant", "content": "..."}`. If None,
            a new conversation is started with a system and user prompt based on
            the `molecule` and selected prompts. Defaults to None.

    Returns:
        tuple[int, str]:
            - `status_code` (int): 200 for success, 400 for failure after retries.
            - `response_text` (str): The raw response text from the LLM if successful,
              otherwise an empty string.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    log_message(f"Calling {LLM} with molecule: {molecule}", logger)

    if detect_seven_member_rings(molecule):
        log_message(f"Detected seven member ring in molecule: {molecule}",
                    logger)
        add_on = ADDON_PROMPT_7_MEMBER
    else:
        add_on = ""

    sys_prompt_final, user_prompt_final, max_completion_tokens = obtain_prompt(
        LLM)
    LLM = LLM.split(":")[0]

    params = {
        "model": LLM,
        "max_completion_tokens": max_completion_tokens,
        "temperature": temperature,
        "seed": 42,
        "top_p": 0.9,
        "metadata": metadata,
    }

    if LLM in DEEPSEEK_MODELS:
        user_prompt_final += add_on

    if "3-7" in LLM:
        params["max_tokens"] = 13192 + 5000
        params["temperature"] = 1
        params.pop("top_p", None)
        params.pop("max_completion_tokens", None)
        params['thinking'] = {"type": "enabled", "budget_tokens": 5000}

    if messages is None:
        messages = [{
            "role": "system",
            "content": sys_prompt_final + add_on
        }, {
            "role":
            "user",
            "content":
            user_prompt_final.replace('{target_smiles}', molecule) + "\n\n" +
            add_on
        }]
    params["messages"] = messages

    try:
        # Call the LLM model
        response = completion(**params)

        res_text = response.choices[0].message.content
    except Exception as e:
        log_message(f"Error in calling {LLM}: {e}", logger)
        log_message(f"Retrying call to {LLM}", logger)
        try:
            response = completion(**params)
            res_text = response.choices[0].message.content
        except Exception as e:
            log_message(f"2nd Error in calling {LLM}: {e}", logger)
            log_message(f"Exiting call to {LLM}", logger)
            return 400, ""
    log_message(f"Received response from LLM: {res_text}", logger)
    return 200, res_text


def split_cot_json(res_text: str) -> tuple[int, list[str], str]:
    """Parses LLM response text to extract Chain-of-Thought (COT) and JSON content.

    Assumes the input `res_text` contains distinct sections demarcated by
    `<cot>\n...thinking steps...</cot>` and `<json>\n...json data...</json>` tags.
    Individual thinking steps within the COT block are expected to be wrapped in
    `<thinking>...</thinking>` tags (though the outer `<thinking>` tag itself seems
    to be part of the split logic rather than the content).

    Args:
        res_text (str):
            The raw response text from an LLM.

    Returns:
        tuple[int, list[str], str]:
            - `status_code` (int): 200 for successful parsing. 501 if COT or JSON
              sections/tags are missing or malformed.
            - `thinking_steps` (list[str]): A list of strings, where each string is
              the content extracted from a `<thinking>...</thinking>` block within the COT.
              Returns an empty list if parsing fails.
            - `json_content` (str): The raw string content extracted from between
              the `<json>\n` and `</json>` tags. Returns an empty string if parsing fails.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        # extract the content within <cot> </cot> tags as thinking content
        thinking_content = res_text[res_text.find("<cot>\n") +
                                    6:res_text.find("</cot>")]
        if not thinking_content:
            return 501, [], ""

        # split the thinking content into individual steps based on the <thinking> </thinking> tags
        thinking_steps = thinking_content.split("<thinking")[1:]
        thinking_steps = [
            step[:step.find("</thinking>")] for step in thinking_steps
        ]
        if not thinking_steps:
            return 501, [], ""
    except Exception as e:
        log_message(f"Error in parsing obtaining COT: {e}", logger)
        return 501, [], ""

    try:
        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
        if not json_content:
            return 501, [], ""
    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 501, [], ""
    return 200, thinking_steps, json_content


def split_json_openAI(res_text: str) -> tuple[int, str]:
    """Extracts JSON content from LLM response text, specifically for OpenAI models.

    This function is designed to parse responses from OpenAI models that are
    expected to provide JSON content, potentially without explicit Chain-of-Thought
    XML-like tags surrounding the JSON block itself. It looks for `<json>\n`
    and `</json>` tags to delimit the JSON content, which implies that even
    OpenAI responses are expected to follow this specific tagging for the JSON part.

    Args:
        res_text (str):
            The raw response text from an OpenAI LLM.

    Returns:
        tuple[int, str]:
            - `status_code` (int): 200 for successful extraction of content between
              `<json>\n` and `</json>` tags. 501 if these tags are not found or
              the content is missing.
            - `json_content` (str): The raw string content extracted from between
              the JSON tags. Returns an empty string if parsing fails.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
        if not json_content:
            return 502, ""

    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 502, ""
    return 200, json_content


def split_json_deepseek(res_text: str) -> tuple[int, list[str], str]:
    """Parses LLM response text from DeepSeek models.

    Assumes DeepSeek response contains a thinking block demarcated by
    `<think>\n...thinking content...</think>` and a JSON block by
    `<json>\n...json data...</json>` tags.

    Args:
        res_text (str):
            The raw response text from a DeepSeek LLM.

    Returns:
        tuple[int, list[str], str]:
            - `status_code` (int): 200 for successful parsing. 503 if tags are
              missing or content is malformed.
            - `thinking_steps` (list[str]): A list containing a single string:
              the content extracted from the `<think>...</think>` block.
              Returns an empty list if parsing fails.
            - `json_content` (str): The raw string content extracted from between
              the `<json>\n` and `</json>` tags. Returns an empty string if parsing fails.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None

    try:
        # extract the content within <cot> </cot> tags as thinking content
        thinking_content = res_text[res_text.find("<think>\n") +
                                    6:res_text.find("</think>")]
        if not thinking_content:
            return 503, [], ""

        json_content = res_text[res_text.find("<json>\n") +
                                7:res_text.find("</json>")]
        if not json_content:
            return 503, [], ""

    except Exception as e:
        log_message(f"Error in parsing LLM response: {e}", logger)
        return 503, [], ""
    return 200, [thinking_content], json_content


def split_json_master(res_text: str, model: str) -> tuple[int, list[str], str]:
    """Master dispatch function for parsing LLM response text based on the model type.

    Selects the appropriate parsing function (`split_json_deepseek`,
    `split_json_openAI`, or `split_cot_json`) based on whether the `model`
    is in `DEEPSEEK_MODELS`, `OPENAI_MODELS`, or defaults to the COT parser.

    Args:
        res_text (str):
            The raw response text from the LLM.
        model (str):
            The identifier of the LLM model used to generate the response (e.g.,
            "deepseek-coder", "gpt-3.5-turbo", "claude-3-opus-20240229").

    Returns:
        tuple[int, list[str], str]:
            - `status_code` (int): Status from the called parsing function, or 505
              if an unexpected error occurs during dispatch.
            - `thinking_steps` (list[str]): List of thinking/COT steps. Empty if
              not applicable (e.g., for OpenAI) or if parsing failed.
            - `json_content` (str): Extracted JSON string. Empty if parsing failed.
    """
    try:
        if model in DEEPSEEK_MODELS:
            status_code, thinking_steps, json_content = split_json_deepseek(
                res_text)
        elif model in OPENAI_MODELS:
            status_code, json_content = split_json_openAI(res_text)
            thinking_steps = []
        else:
            status_code, thinking_steps, json_content = split_cot_json(
                res_text)
    except Exception as e:
        return 505, [], ""

    return status_code, thinking_steps, json_content


def validate_split_json(
        json_content: str) -> tuple[int, list[list[str]] | list[str], list[str], list[float]]:
    """Validates and extracts data from a JSON string (from LLM response).

    Parses the `json_content` string (expected to be valid JSON representation
    of a dictionary) using `ast.literal_eval`. It then attempts to extract
    'data' (pathways/molecules), 'explanation', and 'confidence_scores' keys.

    Args:
        json_content (str):
            The JSON string content extracted from an LLM response.

    Returns:
        tuple[int, list[list[str]] | list[str], list[str], list[float]]:
            - `status_code` (int): 200 for successful parsing and key extraction.
              504 if parsing fails (e.g., invalid JSON) or if expected keys
              are missing.
            - `res_molecules` (list[list[str]] | list[str]): The 'data' field from
              the JSON, typically a list of pathways, where each pathway is a
              list of SMILES strings. Can also be a flat list of SMILES strings
              if the LLM proposes single molecules. Empty list on failure.
            - `res_explanations` (list[str]): The 'explanation' field. Empty list on failure.
            - `res_confidence` (list[float]): The 'confidence_scores' field.
              Empty list on failure.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        result_list = ast.literal_eval(json_content)
        res_molecules = result_list['data']
        res_explanations = result_list['explanation']
        res_confidence = result_list['confidence_scores']
    except Exception as e:
        log_message(f"Error in parsing response: {e}", logger)
        return 504, [], [], []
    return 200, res_molecules, res_explanations, res_confidence


def llm_pipeline(
    molecule: str,
    LLM: str = "claude-3-opus-20240229",
    messages: Optional[list[dict]] = None,
    stability_flag: str = "False",
    hallucination_check: str = "False"
) -> tuple[list[list[str]], list[str], list[float]]:
    """Main pipeline for LLM-based retrosynthesis predictions.

    This function orchestrates the process of:
    1. Calling an LLM (potentially with retries and increasing temperature via `call_LLM`).
       It may switch models (e.g., from DeepSeek to Claude) upon retry.
    2. Parsing the LLM's response using `split_json_master`.
    3. Validating the extracted JSON content using `validate_split_json`.
    4. (Implicitly, further down in the full code) Applying validity checks,
       stability checks (`stability_checker`), and hallucination checks
       (`hallucination_checker`) to the results.

    The pipeline attempts to get valid pathways. The number of retries is
    controlled by `max_run` which is influenced by `stability_flag` and
    `hallucination_check`.

    Args:
        molecule (str):
            The SMILES string of the target molecule.
        LLM (str, optional):
            The primary LLM identifier. Defaults to "claude-3-opus-20240229".
        messages (Optional[list[dict]], optional):
            Conversation history for the LLM. Defaults to None (new conversation).
        stability_flag (str, optional):
            "True" or "False" to enable stability checks. Affects retry attempts.
            Defaults to "False".
        hallucination_check (str, optional):
            "True" or "False" to enable hallucination checks. Affects retry attempts.
            Defaults to "False".

    Returns:
        tuple[list[list[str]], list[str], list[float]]:
            - `output_pathways` (list[list[str]]): A list of proposed synthetic
              pathways. Each pathway is a list of SMILES strings (precursors).
              Returns an empty list if no valid pathways are found after all attempts.
            - `output_explanations` (list[str]): Corresponding explanations for each pathway.
            - `output_confidence` (list[float]): Corresponding confidence scores for each pathway.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    output_pathways: list[list[str]] = []
    output_explanations: list[str] = []
    output_confidence: list[float] = []
    run = 0.0
    if stability_flag.lower() == "true" or hallucination_check.lower(
    ) == "true":
        max_run = 1.5
    else:
        max_run = 0.6
    while (output_pathways == [] and run < max_run):
        log_message(f"Calling LLM with molecule: {molecule} and run: {run}",
                    logger)

        # Selecting the model based on the run number
        current_model = LLM
        if LLM in DEEPSEEK_MODELS and run > 0.0:
            current_model = "claude-3-opus-20240229"

        # --------------------
        # Call LLM
        status_code, res_text = call_LLM(molecule,
                                         current_model,
                                         messages=messages,
                                         temperature=run)
        if status_code != 200:
            log_message(f"Error in calling LLM: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Split the response text
        status_code, thinking_steps, json_content = split_json_master(
            res_text, current_model)
        if status_code != 200:
            log_message(f"Error in splitting cot json: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Validate the split json content
        status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
            json_content)
        if status_code != 200:
            log_message(f"Error in validating split json content: {res_text}",
                        logger)
            run += 0.1
            get_error_log(status_code)
            continue

        # --------------------
        # Check the validity of the molecules obtained from LLM
        output_pathways, output_explanations, output_confidence = validity_check(
            molecule, res_molecules, res_explanations, res_confidence)

        # --------------------
        # Stability check
        if stability_flag.lower() == "true":
            status_code, stable_pathways = stability_checker(output_pathways)
            if status_code != 200:
                log_message(f"Error in stability check: {stable_pathways}",
                            logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = stable_pathways

        # --------------------
        # Hallucination check
        if hallucination_check.lower() == "true":
            log_message(f"Calling hallucination check with pathways: {output_pathways}",
                        logger)
            status_code, hallucination_pathways = hallucination_checker(
                molecule, output_pathways)
            if status_code != 200:
                log_message(
                    f"Error in hallucination check: {hallucination_pathways}",
                    logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = hallucination_pathways

        # --------------------
        # Update the messages for the next call
        # if output_pathways == []:
        #     messages = [{
        #         "role": "system",
        #         "content": SYS_PROMPT
        #     }, {
        #         "role": "user",
        #         "content": USER_PROMPT
        #     }, {
        #         "role": "assistant",
        #         "content": res_text
        #     }, {
        #         "role": "user",
        #         "content": "<add something here>"
        #     }]

        log_message(
            f"Output Pathways: {output_pathways},\n\
                Output Explanations: {output_explanations},\n\
                    Output Confidence: {output_confidence}", logger)
        run += 0.1

    return output_pathways, output_explanations, output_confidence


def get_error_log(status_code: int) -> None:
    """Logs a descriptive error message based on a status code.

    Retrieves an error description from `src.variables.ERROR_MAP` using the
    provided `status_code` and logs it. If the status code is not found
    in the map, it logs an "unrecognized" error message.

    Args:
        status_code (int):
            The status code representing an error that occurred.
    """
    logger = context_logger.get() if ENABLE_LOGGING else None

    if status_code in ERROR_MAP:
        description = ERROR_MAP[status_code]
        log_message(f"Error Code: {status_code},\n Description: {description}",
                    logger)
    else:
        log_message(f"Error Code: {status_code} is not recognized.", logger)
