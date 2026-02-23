"""LLM-based retrosynthesis: prompts, completion, and response parsing.

Calls LLMs (Claude, GPT, DeepSeek via LiteLLM) to predict retrosynthetic
precursors. Parses chain-of-thought and JSON responses, validates outputs,
and runs optional stability/hallucination checks.
"""

import ast
import os
import re
from typing import Optional

import litellm
from dotenv import load_dotenv
from litellm import completion

from src.cache import cache_results
from src.protecting_group import mask_protecting_groups_multisymbol
from src.utils.hallucination_checks import hallucination_checker
from src.utils.job_context import logger as context_logger
from src.utils.langfuse_config import get_langfuse_metadata
from src.utils.stability_checks import stability_checker
from src.utils.utils_molecule import detect_seven_member_rings, validity_check
from src.variables import (
    ADDON_PROMPT_7_MEMBER,
    DEEPSEEK_MODELS,
    ERROR_MAP,
    OPENAI_MODELS,
    PROTECTING_GROUP_CONTEXT,
    SYS_PROMPT,
    SYS_PROMPT_DEEPSEEK,
    SYS_PROMPT_OPENAI,
    SYS_PROMPT_V4,
    USER_PROMPT,
    USER_PROMPT_DEEPSEEK,
    USER_PROMPT_DEEPSEEK_V4,
    USER_PROMPT_OPENAI,
    USER_PROMPT_V4,
)

load_dotenv()

litellm.success_callback = ["langfuse"]
litellm.drop_params = True

ENABLE_LOGGING = (
    True if os.getenv("ENABLE_LOGGING", "false").lower() == "true" else False
)

_MAX_API_RETRIES = 2  # Number of API attempts before returning 400
_TEMPERATURE_STEP = 0.1  # Temperature increment per retry attempt
_FALLBACK_MODEL = "claude-opus-4-6"  # Fallback when DeepSeek fails on retry

# (model_family, advanced) -> (sys_prompt, user_prompt, max_completion_tokens)
_PROMPT_CONFIG = {
    ("deepseek", False): (SYS_PROMPT_DEEPSEEK, USER_PROMPT_DEEPSEEK, 16384),
    ("deepseek", True): (SYS_PROMPT_V4, USER_PROMPT_DEEPSEEK_V4, 16384),
    ("openai", False): (SYS_PROMPT_OPENAI, USER_PROMPT_OPENAI, 8192),
    ("openai", True): (SYS_PROMPT_OPENAI, USER_PROMPT_OPENAI, 8192),
    ("default", False): (SYS_PROMPT, USER_PROMPT, 4096),
    ("default", True): (SYS_PROMPT_V4, USER_PROMPT_V4, 4096),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_logger():
    """Return the context logger if logging is enabled, else None."""
    if ENABLE_LOGGING:
        return context_logger.get(None)
    return None


def _log(message: str, logger=None):
    """Log via the provided logger, falling back to print.

    Parameters
    ----------
    message : str
        The message to be logged
    logger : Logger, optional
        The logger object, by default None

    Returns
    -------
    None
    """
    if ENABLE_LOGGING:
        if logger is not None:
            logger.info(message)
        else:
            print(message)
    else:
        pass


def _log_error(status_code: int):
    """Log a descriptive error message for the given status code.

    Parameters
    ----------
    status_code : int
        The status code to log

    Returns
    -------
    None
    """
    if ENABLE_LOGGING:
        logger = _get_logger()
        description = ERROR_MAP.get(status_code, "Unrecognized error")
        _log(f"Error {status_code}: {description}", logger)


def extract_tag_content(text: str, tag: str) -> str | None:
    """Extract content between ``<tag>`` and ``</tag>``. Returns *None* if not found.

    Parameters
    ----------
    text : str
        The text to extract the tag from
    tag : str
        The tag to extract

    Returns
    -------
    str | None

    Examples
    --------
    >>> extract_tag_content("<json>hello</json>", "json")
    'hello'
    >>> extract_tag_content("<tag id='x'>content</tag>", "tag")
    'content'
    >>> extract_tag_content("no tags here", "json") is None
    True
    """
    match = re.search(
        rf"<{re.escape(tag)}\b[^>]*>\s*(.*?)\s*</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def _classify_model(model: str) -> str:
    """Classify a model identifier into 'deepseek', 'openai', or 'default'.

    Parameters
    ----------
    model : str
        The model identifier to classify

    Returns
    -------
    str
    """
    if model in DEEPSEEK_MODELS:
        return "deepseek"
    if model in OPENAI_MODELS:
        return "openai"
    return "default"


# ---------------------------------------------------------------------------
# Prompt selection
# ---------------------------------------------------------------------------


def obtain_prompt(model: str) -> tuple[str, str, int]:
    """Select system prompt, user prompt, and max tokens for *model*.

    Append ``:adv`` to the model string (e.g. ``claude-opus-4-6:adv``) to
    select the advanced prompt variant.

    Parameters
    ----------
    model : str
        The model identifier to classify

    Returns
    -------
    tuple[str, str, int]
        The system prompt, user prompt, and max tokens

    Examples
    --------
    >>> sys_prompt, user_prompt, max_tokens = obtain_prompt("gpt-4o")
    >>> max_tokens
    8192
    >>> len(sys_prompt) > 0 and len(user_prompt) > 0
    True
    >>> _, _, tokens = obtain_prompt("claude-opus-4-6:adv")
    >>> tokens
    4096
    """
    parts = model.split(":")
    advanced = len(parts) > 1 and parts[1] == "adv"
    family = _classify_model(parts[0])
    return _PROMPT_CONFIG[(family, advanced)]


# ---------------------------------------------------------------------------
# LLM call
# ---------------------------------------------------------------------------


def build_addon_prompt(
    molecule: str,
    use_protecting_group_feature: bool,
    logger,
) -> str:
    """Build supplementary prompt text for seven-member rings / protecting groups.

    Appends extra instructions when the molecule contains seven-member rings
    or when protecting-group masking is enabled.

    Parameters
    ----------
    molecule : str
        The molecule to build the addon prompt for
    use_protecting_group_feature : bool
        Whether to use the protecting group feature
    logger : Logger
        The logger to use

    Returns
    -------
    str
        The addon prompt
    """
    parts: list[str] = []

    if detect_seven_member_rings(molecule):
        _log(f"Detected seven-member ring in molecule: {molecule}", logger)
        parts.append(ADDON_PROMPT_7_MEMBER)

    if use_protecting_group_feature:
        masked = mask_protecting_groups_multisymbol(molecule)
        if masked not in (molecule, "INVALID_SMILES"):
            _log(f"Detected protecting groups: {molecule} -> {masked}", logger)
            parts.append(
                PROTECTING_GROUP_CONTEXT.format(molecule=molecule, masked_smiles=masked)
            )

    return "".join(parts)


def build_completion_params(
    model: str,
    messages: list[dict],
    max_completion_tokens: int,
    temperature: float,
) -> dict:
    """Assemble the kwargs dict for ``litellm.completion``.

    Adds seed, and Langfuse metadata. Special handling for
    extended-thinking models (e.g. ``3-7`` in model name).

    Parameters
    ----------
    model : str
        The model to build the completion parameters for
    messages : list[dict]
        The messages to build the completion parameters for
    max_completion_tokens : int
        The max completion tokens to build the completion parameters for
    temperature : float
        The temperature to build the completion parameters for

    Returns
    -------
    dict
        The completion parameters
    """
    params: dict = {
        "model": model,
        "messages": messages,
        "max_completion_tokens": max_completion_tokens,
        "temperature": temperature,
        "seed": 42,
        "metadata": get_langfuse_metadata("retrosynthesis"),
    }

    if "3-7" in model:
        params["max_tokens"] = 18192
        params["temperature"] = 1
        params.pop("max_completion_tokens", None)
        params["thinking"] = {"type": "enabled", "budget_tokens": 5000}

    return params


@cache_results
def call_LLM(
    molecule: str,
    model: str = "claude-opus-4-6",
    temperature: float = 0.0,
    messages: Optional[list[dict]] = None,
    use_protecting_group_feature: bool = False,
) -> tuple[int, str]:
    """Call an LLM to predict retrosynthetic precursors for *molecule*.

    Returns 200 on success, 400 on failure.

    Example
    -------
    >>> from deepretro.utils.llm import call_LLM
    >>> from src.cache import clear_cache_for_molecule
    >>> clear_cache_for_molecule("C1CCCCC1")
    >>> status, res_text = call_LLM("C1CCCCC1", model="claude-haiku-4-5-20251001")
    >>> print(status)
    200

    Parameters
    ----------
    molecule : str
        The molecule to call the LLM for
    model : str
        The model to call the LLM for
    temperature : float
        The temperature to call the LLM for
    messages : list[dict], optional
        The messages to call the LLM for
    use_protecting_group_feature : bool, optional
        Whether to use the protecting group feature

    Returns
    -------
    tuple[int, str]
        The status code and the response text
    """
    logger = _get_logger()
    _log(f"Calling {model} with molecule: {molecule}", logger)

    addon = build_addon_prompt(molecule, use_protecting_group_feature, logger)
    sys_prompt, user_prompt, max_tokens = obtain_prompt(model)
    model_name = model.split(":")[0]

    if model_name in DEEPSEEK_MODELS:
        user_prompt += addon

    if messages is None:
        messages = [
            {"role": "system", "content": sys_prompt + addon},
            {
                "role": "user",
                "content": (
                    user_prompt.replace("{target_smiles}", molecule) + "\n\n" + addon
                ),
            },
        ]

    params = build_completion_params(model_name, messages, max_tokens, temperature)

    for attempt in range(1, _MAX_API_RETRIES + 1):
        try:
            response = completion(**params)
            res_text = response.choices[0].message.content
            _log(f"Received response from LLM: {res_text}", logger)
            return 200, res_text
        except Exception as exc:
            _log(
                f"Attempt {attempt}/{_MAX_API_RETRIES} failed for {model_name}: {exc}",
                logger,
            )

    _log(f"All {_MAX_API_RETRIES} attempts failed for {model_name}", logger)
    return 400, ""


# ---------------------------------------------------------------------------
# Response parsing
# ---------------------------------------------------------------------------


def _parse_cot(res_text: str) -> tuple[int, list[str], str]:
    """Parse a chain-of-thought response (``<cot>`` / ``<thinking>`` / ``<json>``).

    Parameters
    ----------
    res_text : str
        The response text to parse

    Returns
    -------
    tuple[int, list[str], str]
        The status code, the thinking steps, and the JSON content
    """
    cot_content = extract_tag_content(res_text, "cot")
    if not cot_content:
        return 501, [], ""

    thinking_steps = [
        m.group(1)
        for m in re.finditer(r"<thinking[^>]*>(.*?)</thinking>", cot_content, re.DOTALL)
    ]
    if not thinking_steps:
        return 501, [], ""

    json_content = extract_tag_content(res_text, "json")
    if not json_content:
        return 501, [], ""

    return 200, thinking_steps, json_content


def _parse_deepseek(res_text: str) -> tuple[int, list[str], str]:
    """Parse a DeepSeek response (``<think>`` / ``<json>``).

    Parameters
    ----------
    res_text : str
        The response text to parse

    Returns
    -------
    tuple[int, list[str], str]
        The status code, the thinking content, and the JSON content
    """
    thinking_content = extract_tag_content(res_text, "think")
    if not thinking_content:
        return 503, [], ""

    json_content = extract_tag_content(res_text, "json")
    if not json_content:
        return 503, [], ""

    return 200, [thinking_content], json_content


def parse_response(res_text: str, model: str) -> tuple[int, list[str], str]:
    """Parse an LLM response into ``(status_code, thinking_steps, json_content)``.

    Dispatches to the correct strategy based on model family.

    Parameters
    ----------
    res_text : str
        The response text to parse
    model : str
        The model to parse the response for

    Returns
    -------
    tuple[int, list[str], str]
        The status code, the thinking steps, and the JSON content

    Examples
    --------
    >>> res = '<cot><thinking>Step 1</thinking></cot><json>{"data":[["CC"]],"explanation":["x"],"confidence_scores":[0.9]}</json>'
    >>> status, steps, json_str = parse_response(res, "claude-opus-4-6")
    >>> status
    200
    >>> steps
    ['Step 1']
    >>> 'CC' in json_str
    True
    >>> res_openai = '<json>{"data":[["CC"]],"explanation":["x"],"confidence_scores":[0.9]}</json>'
    >>> status, steps, _ = parse_response(res_openai, "gpt-4o")
    >>> status
    200
    >>> steps
    []
    >>> parse_response("", "claude-opus-4-6")[0]
    501
    """
    logger = _get_logger()
    family = _classify_model(model)

    try:
        if family == "deepseek":
            return _parse_deepseek(res_text)

        if family == "openai":
            json_content = extract_tag_content(res_text, "json")
            if not json_content:
                return 502, [], ""
            return 200, [], json_content

        return _parse_cot(res_text)
    except Exception as exc:
        _log(f"Error parsing {family} response: {exc}", logger)
        return 505, [], ""


def validate_json_response(
    json_content: str,
) -> tuple[int, list[str], list[str], list[int]]:
    """Parse JSON content and extract molecules, explanations, and confidence scores.

    Parameters
    ----------
    json_content : str
        The JSON content to validate

    Returns
    -------
    tuple[int, list[str], list[str], list[int]]
        The status code, the molecules, the explanations, and the confidence scores

    Examples
    --------
    >>> json_str = '{"data":[["CC(=O)C","CC=O"],["CC(O)CC"]],"explanation":["Aldol","Oxidation"],"confidence_scores":[0.9,0.8]}'
    >>> status, data, expl, conf = validate_json_response(json_str)
    >>> status
    200
    >>> data
    [['CC(=O)C', 'CC=O'], ['CC(O)CC']]
    >>> expl
    ['Aldol', 'Oxidation']
    >>> conf
    [0.9, 0.8]
    """
    logger = _get_logger()
    try:
        result = ast.literal_eval(json_content)
        return 200, result["data"], result["explanation"], result["confidence_scores"]
    except Exception as exc:
        _log(f"Error parsing JSON response: {exc}", logger)
        return 504, [], [], []


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def llm_pipeline(
    molecule: str,
    model: str = "claude-opus-4-6",
    messages: Optional[list[dict]] = None,
    stability_check: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
) -> tuple[list[list[str]], list[str], list[float]]:
    """End-to-end retrosynthesis pipeline.

    Calls the LLM, parses and validates the response, and optionally runs
    stability / hallucination checks.  Retries with increasing temperature on
    failure, falling back to Claude when DeepSeek models fail on retry.

    Example
    -------
    >>> from deepretro.utils.llm import llm_pipeline
    >>> pathways, explanations, confidence = llm_pipeline("C1CCCCC1")  # doctest: +SKIP
    >>> isinstance(pathways, list) and isinstance(explanations, list)  # doctest: +SKIP
    True

    Parameters
    ----------
    molecule : str
        The molecule to run the pipeline for
    model : str
        The model to run the pipeline for
    messages : list[dict], optional
        The messages to run the pipeline for
    stability_check : str, optional
        Whether to run the stability check
    hallucination_check : str, optional
        Whether to run the hallucination check
    use_protecting_group_feature : bool, optional
        Whether to use the protecting group feature

    Returns
    -------
    tuple[list[list[str]], list[str], list[float]]
        The output pathways, explanations, and confidence scores
    """
    logger = _get_logger()
    max_attempts = (
        15
        if (stability_check.lower() == "true" or hallucination_check.lower() == "true")
        else 6
    )

    output_pathways: list[list[str]] = []
    output_explanations: list[str] = []
    output_confidence: list[float] = []

    for attempt in range(max_attempts):
        temperature = round(attempt * _TEMPERATURE_STEP, 2)
        current_model = model
        if model in DEEPSEEK_MODELS and attempt > 0:
            current_model = _FALLBACK_MODEL

        status, res_text = call_LLM(
            molecule,
            current_model,
            messages=messages,
            temperature=temperature,
            use_protecting_group_feature=use_protecting_group_feature,
        )
        if status != 200:
            _log_error(status)
            continue

        status, _thinking_steps, json_content = parse_response(res_text, current_model)
        if status != 200:
            _log_error(status)
            continue

        status, res_molecules, res_explanations, res_confidence = (
            validate_json_response(json_content)
        )
        if status != 200:
            _log_error(status)
            continue

        output_pathways, output_explanations, output_confidence = validity_check(
            molecule, res_molecules, res_explanations, res_confidence
        )

        if stability_check.lower() == "true" and output_pathways:
            status, stable_pathways = stability_checker(output_pathways)
            if status != 200:
                _log_error(status)
                output_pathways = []
                continue
            output_pathways = stable_pathways

        if hallucination_check.lower() == "true" and output_pathways:
            _log(
                f"Running hallucination check on {len(output_pathways)} pathways",
                logger,
            )
            status, checked_pathways = hallucination_checker(molecule, output_pathways)
            if status != 200:
                _log_error(status)
                output_pathways = []
                continue
            output_pathways = checked_pathways

        if output_pathways:
            _log(
                f"Success: {len(output_pathways)} pathways, "
                f"{len(output_explanations)} explanations",
                logger,
            )
            break

    return output_pathways, output_explanations, output_confidence
