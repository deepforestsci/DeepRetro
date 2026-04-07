"""LLM retrosynthesis pipeline (ported from src.utils.llm)."""

import os
import ast
import litellm
from typing import Optional
from dotenv import load_dotenv
from litellm import completion
from deepretro.utils.variables import OPENAI_MODELS, DEEPSEEK_MODELS
from deepretro.utils.variables import USER_PROMPT, SYS_PROMPT
from deepretro.utils.variables import USER_PROMPT_V4, SYS_PROMPT_V4
from deepretro.utils.variables import USER_PROMPT_OPENAI, SYS_PROMPT_OPENAI
from deepretro.utils.variables import USER_PROMPT_DEEPSEEK, SYS_PROMPT_DEEPSEEK
from deepretro.utils.variables import ADDON_PROMPT_7_MEMBER, USER_PROMPT_DEEPSEEK_V4
from deepretro.utils.variables import ERROR_MAP, PROTECTING_GROUP_CONTEXT
from deepretro.utils.utils_molecule import validity_check, detect_seven_member_rings
from deepretro.migration.job_context import logger as context_logger
from deepretro.migration.pipeline_checks import stability_checker, hallucination_checker
from deepretro.migration.protecting_group import mask_protecting_groups_multisymbol

load_dotenv()

litellm.success_callback = ["langfuse"]
litellm.drop_params = True

from deepretro.migration.langfuse_config import get_langfuse_metadata

ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING",
                                    "true").lower() == "false" else True


def log_message(message: str, logger=None):
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def obtain_prompt(LLM: str):
    """Obtain the prompt based on the LLM model.

    Returns
    -------
    str, str, int
        The system prompt, user prompt and max completion tokens.
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


def call_LLM(molecule: str,
             LLM: str = "claude-opus-4-20250514",
             temperature: float = 0.0,
             messages: Optional[list[dict]] = None,
             use_protecting_group_feature: bool = False) -> tuple[int, str]:
    """Call the LLM model to predict the next retrosynthetic step."""
    logger = context_logger.get() if ENABLE_LOGGING else None
    log_message(f"Calling {LLM} with molecule: {molecule}", logger)

    if detect_seven_member_rings(molecule):
        log_message(f"Detected seven member ring in molecule: {molecule}",
                    logger)
        add_on = ADDON_PROMPT_7_MEMBER
    else:
        add_on = ""

    if use_protecting_group_feature:
        masked_smiles = mask_protecting_groups_multisymbol(molecule)
        if masked_smiles != molecule and masked_smiles != "INVALID_SMILES":
            log_message(
                f"Detected protecting groups in molecule: {molecule} -> {masked_smiles}",
                logger)
            protecting_group_context = PROTECTING_GROUP_CONTEXT.format(
                molecule=molecule, masked_smiles=masked_smiles)
            add_on += protecting_group_context

    sys_prompt_final, user_prompt_final, max_completion_tokens = obtain_prompt(
        LLM)
    is_thinking = "think" in LLM.split(":")
    LLM = LLM.split(":")[0]

    params = {
        "model": LLM,
        "max_completion_tokens": max_completion_tokens,
        "temperature": temperature,
        "seed": 42,
        "top_p": 0.9,
        "metadata": get_langfuse_metadata("retrosynthesis"),
    }

    if LLM in DEEPSEEK_MODELS:
        user_prompt_final += add_on

    if is_thinking:
        params["max_tokens"] = 16384
        params.pop("temperature", None)
        params.pop("top_p", None)
        params.pop("seed", None)
        params.pop("max_completion_tokens", None)
        params['thinking'] = {"type": "adaptive"}

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
    """Parse LLM response to extract thinking steps and json content."""
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        thinking_content = res_text[res_text.find("<cot>\n") +
                                    6:res_text.find("</cot>")]
        if not thinking_content:
            return 501, [], ""

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
    """Split response from OpenAI models (no COT)."""
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
    """Parse DeepSeek response to extract thinking and json content."""
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
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
    """Route response parsing based on model type."""
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
    except Exception:
        return 505, [], ""

    return status_code, thinking_steps, json_content


def validate_split_json(
        json_content: str) -> tuple[int, list[str], list[str], list[int]]:
    """Validate parsed json content from LLM response."""
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
    LLM: str = "claude-opus-4-20250514",
    messages: Optional[list[dict]] = None,
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
    hallucination_checker_fn=None,
) -> tuple[list[list[str]], list[str], list[float]]:
    """Full pipeline: call LLM, parse, validate, stability/hallucination check."""
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

        current_model = LLM
        if LLM in DEEPSEEK_MODELS and run > 0.0:
            current_model = "claude-opus-4-20250514"

        status_code, res_text = call_LLM(
            molecule,
            current_model,
            messages=messages,
            temperature=run,
            use_protecting_group_feature=use_protecting_group_feature)
        if status_code != 200:
            log_message(f"Error in calling LLM: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        status_code, thinking_steps, json_content = split_json_master(
            res_text, current_model)
        if status_code != 200:
            log_message(f"Error in splitting cot json: {res_text}", logger)
            run += 0.1
            get_error_log(status_code)
            continue

        status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
            json_content)
        if status_code != 200:
            log_message(f"Error in validating split json content: {res_text}",
                        logger)
            run += 0.1
            get_error_log(status_code)
            continue

        output_pathways, output_explanations, output_confidence = validity_check(
            molecule, res_molecules, res_explanations, res_confidence)

        if stability_flag.lower() == "true":
            status_code, stable_pathways = stability_checker(output_pathways)
            if status_code != 200:
                log_message(f"Error in stability check: {stable_pathways}",
                            logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = stable_pathways

        if hallucination_check.lower() == "true":
            log_message(
                f"Calling hallucination check with pathways: {output_pathways}",
                logger)
            _checker = hallucination_checker_fn or hallucination_checker
            status_code, hallucination_pathways = _checker(
                molecule, output_pathways)
            if status_code != 200:
                log_message(
                    f"Error in hallucination check: {hallucination_pathways}",
                    logger)
                run += 0.1
                get_error_log(status_code)
                continue
            output_pathways = hallucination_pathways

        log_message(
            f"Output Pathways: {output_pathways},\n\
                Output Explanations: {output_explanations},\n\
                    Output Confidence: {output_confidence}", logger)
        run += 0.1

    return output_pathways, output_explanations, output_confidence


def get_error_log(status_code: int):
    """Print error message based on status code."""
    logger = context_logger.get() if ENABLE_LOGGING else None

    if status_code in ERROR_MAP:
        description = ERROR_MAP[status_code]
        log_message(f"Error Code: {status_code},\n Description: {description}",
                    logger)
    else:
        log_message(f"Error Code: {status_code} is not recognized.", logger)
