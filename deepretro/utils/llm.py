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

ENABLE_LOGGING = os.getenv("ENABLE_LOGGING", "true").lower() != "false"

_MAX_API_RETRIES = 2
_TEMPERATURE_STEP = 0.1
_FALLBACK_MODEL = "claude-opus-4-6"

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
    return context_logger.get() if ENABLE_LOGGING else None


def _log(message: str, logger=None):
    """Log via the provided logger, falling back to print."""
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def _log_error(status_code: int):
    """Log a descriptive error message for the given status code."""
    logger = _get_logger()
    description = ERROR_MAP.get(status_code, "Unrecognized error")
    _log(f"Error {status_code}: {description}", logger)


def _extract_tag(text: str, tag: str) -> str | None:
    """Extract content between ``<tag>`` and ``</tag>``. Returns *None* if not found."""
    match = re.search(
        rf"<{re.escape(tag)}\b[^>]*>\s*(.*?)\s*</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def _classify_model(model: str) -> str:
    """Classify a model identifier into 'deepseek', 'openai', or 'default'."""
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
    """
    parts = model.split(":")
    advanced = len(parts) > 1 and parts[1] == "adv"
    family = _classify_model(parts[0])
    return _PROMPT_CONFIG[(family, advanced)]


# ---------------------------------------------------------------------------
# LLM call
# ---------------------------------------------------------------------------


def _build_addon_prompt(
    molecule: str,
    use_protecting_group_feature: bool,
    logger,
) -> str:
    """Build supplementary prompt text for seven-member rings / protecting groups."""
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


def _build_completion_params(
    model: str,
    messages: list[dict],
    max_completion_tokens: int,
    temperature: float,
) -> dict:
    """Assemble the kwargs dict for ``litellm.completion``."""
    params: dict = {
        "model": model,
        "messages": messages,
        "max_completion_tokens": max_completion_tokens,
        "temperature": temperature,
        "seed": 42,
        "top_p": 0.9,
        "metadata": get_langfuse_metadata("retrosynthesis"),
    }

    if "3-7" in model:
        params["max_tokens"] = 18192
        params["temperature"] = 1
        params.pop("top_p", None)
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

    Returns ``(status_code, response_text)`` — 200 on success, 400 on failure.
    """
    logger = _get_logger()
    _log(f"Calling {model} with molecule: {molecule}", logger)

    addon = _build_addon_prompt(molecule, use_protecting_group_feature, logger)
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

    params = _build_completion_params(model_name, messages, max_tokens, temperature)

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
    """Parse a chain-of-thought response (``<cot>`` / ``<thinking>`` / ``<json>``)."""
    cot_content = _extract_tag(res_text, "cot")
    if not cot_content:
        return 501, [], ""

    thinking_steps = [
        m.group(1)
        for m in re.finditer(r"<thinking[^>]*>(.*?)</thinking>", cot_content, re.DOTALL)
    ]
    if not thinking_steps:
        return 501, [], ""

    json_content = _extract_tag(res_text, "json")
    if not json_content:
        return 501, [], ""

    return 200, thinking_steps, json_content


def _parse_deepseek(res_text: str) -> tuple[int, list[str], str]:
    """Parse a DeepSeek response (``<think>`` / ``<json>``)."""
    thinking_content = _extract_tag(res_text, "think")
    if not thinking_content:
        return 503, [], ""

    json_content = _extract_tag(res_text, "json")
    if not json_content:
        return 503, [], ""

    return 200, [thinking_content], json_content


def parse_response(res_text: str, model: str) -> tuple[int, list[str], str]:
    """Parse an LLM response into ``(status_code, thinking_steps, json_content)``.

    Dispatches to the correct strategy based on model family.
    """
    logger = _get_logger()
    family = _classify_model(model)

    try:
        if family == "deepseek":
            return _parse_deepseek(res_text)

        if family == "openai":
            json_content = _extract_tag(res_text, "json")
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
    """Parse JSON content and extract molecules, explanations, and confidence scores."""
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
    stability_check: bool = False,
    hallucination_check: bool = False,
    use_protecting_group_feature: bool = False,
) -> tuple[list[list[str]], list[str], list[float]]:
    """End-to-end retrosynthesis pipeline.

    Calls the LLM, parses and validates the response, and optionally runs
    stability / hallucination checks.  Retries with increasing temperature on
    failure, falling back to Claude when DeepSeek models fail on retry.
    """
    logger = _get_logger()
    max_attempts = 15 if (stability_check or hallucination_check) else 6

    output_pathways: list[list[str]] = []
    output_explanations: list[str] = []
    output_confidence: list[float] = []

    for attempt in range(max_attempts):
        temperature = round(attempt * _TEMPERATURE_STEP, 2)
        current_model = model
        if model in DEEPSEEK_MODELS and attempt > 0:
            current_model = _FALLBACK_MODEL

        _log(
            f"Pipeline attempt {attempt + 1}/{max_attempts}, "
            f"temp={temperature:.1f}, model={current_model}",
            logger,
        )

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

        if stability_check and output_pathways:
            status, stable_pathways = stability_checker(output_pathways)
            if status != 200:
                _log_error(status)
                output_pathways = []
                continue
            output_pathways = stable_pathways

        if hallucination_check and output_pathways:
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
