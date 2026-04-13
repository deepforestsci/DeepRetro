"""LLM-based single-step retrosynthesis helpers.

This module calls chat-completion LLMs through LiteLLM, parses tagged model
responses, and filters candidate precursor pathways. Provider-specific model
normalization and completion argument construction live in
``deepretro.utils.llm_helpers`` so the public pipeline stays focused on the
retrosynthesis workflow.
"""

from __future__ import annotations

import ast
import re
import time
from typing import Any, cast

import litellm
import structlog
from dotenv import load_dotenv
from litellm import completion

from deepretro.algorithms.hallucination_checker import calculate_hallucination_score
from deepretro.algorithms.stability_checker import check_molecule_stability
from deepretro.utils.llm_helpers import (
    ChatMessage,
    Pathway,
    PromptMode,
    ThinkingEffort,
    build_completion_params,
    coerce_response_text,
    extract_json_payload,
    extract_tag_content,
    is_enabled,
    resolve_model_selection,
)
from deepretro.utils.utils_molecule import detect_seven_member_rings, validity_check
from deepretro.utils.variables import (
    ADDON_PROMPT_7_MEMBER,
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

logger = structlog.get_logger(__name__)

MAX_API_RETRIES = 2
TEMPERATURE_STEP = 0.1
FALLBACK_MODEL = "claude-opus-4-6"
STABLE_ASSESSMENTS = {"Likely stable", "Moderately stable"}
HALLUCINATION_SEVERITIES = {"low", "medium"}

PROMPT_CONFIG = {
    ("deepseek", "standard"): (SYS_PROMPT_DEEPSEEK, USER_PROMPT_DEEPSEEK, 16384),
    ("deepseek", "advanced"): (SYS_PROMPT_V4, USER_PROMPT_DEEPSEEK_V4, 16384),
    ("openai", "standard"): (SYS_PROMPT_OPENAI, USER_PROMPT_OPENAI, 16384),
    ("openai", "advanced"): (SYS_PROMPT_OPENAI, USER_PROMPT_OPENAI, 16384),
    ("default", "standard"): (SYS_PROMPT, USER_PROMPT, 16384),
    ("default", "advanced"): (SYS_PROMPT_V4, USER_PROMPT_V4, 16384),
}


def obtain_prompt(
    model: str,
    prompt_mode: PromptMode | None = None,
) -> tuple[str, str, int]:
    """Select the prompt pair and default output-token limit for a model.

    Parameters
    ----------
    model : str
        Model identifier, optionally with ``:adv``.
    prompt_mode : {"standard", "advanced"}, optional
        Explicit prompt-mode override. When omitted, ``:adv`` is honored.

    Returns
    -------
    tuple[str, str, int]
        System prompt, user prompt, and default output-token limit.

    Examples
    --------
    >>> _, _, max_tokens = obtain_prompt("openai/gpt-4o-mini")
    >>> max_tokens
    16384
    """
    selection = resolve_model_selection(model, prompt_mode=prompt_mode)
    return PROMPT_CONFIG[(selection.family, selection.prompt_mode)]


def build_addon_prompt(molecule: str) -> str:
    """Build non-provider prompt context for a target molecule.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.

    Returns
    -------
    str
        Additional prompt text. Currently this only adds seven-membered ring
        guidance.

    Examples
    --------
    >>> build_addon_prompt("CCO")
    ''
    """
    if detect_seven_member_rings(molecule):
        logger.debug("Detected seven-membered ring", molecule=molecule)
        return ADDON_PROMPT_7_MEMBER
    return ""


def build_messages(
    molecule: str,
    model: str,
    addon: str = "",
    messages: list[ChatMessage] | None = None,
    prompt_mode: PromptMode | None = None,
) -> list[ChatMessage]:
    """Return caller-provided messages or build the default prompt pair.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.
    model : str
        LiteLLM model identifier used for prompt selection.
    addon : str, optional
        Additional context appended to the prompt.
    messages : list[ChatMessage], optional
        Caller-provided messages. When provided, these are returned unchanged.
    prompt_mode : {"standard", "advanced"}, optional
        Prompt-mode override.

    Returns
    -------
    list[ChatMessage]
        Chat messages ready for ``litellm.completion``.

    Examples
    --------
    >>> messages = build_messages("CCO", "openai/gpt-4o-mini")
    >>> [message["role"] for message in messages]
    ['system', 'user']
    >>> custom = [{"role": "user", "content": "Reply OK"}]
    >>> build_messages("CCO", "openai/gpt-4o-mini", messages=custom) is custom
    True
    """
    if messages is not None:
        return messages

    sys_prompt, user_prompt, _ = obtain_prompt(model, prompt_mode=prompt_mode)
    user_content = user_prompt.replace("{target_smiles}", molecule)
    if addon:
        user_content = f"{user_content}\n\n{addon}"

    return [
        {"role": "system", "content": sys_prompt + addon},
        {"role": "user", "content": user_content},
    ]


def call_LLM(
    molecule: str,
    model: str = "claude-opus-4-6",
    temperature: float = 0.0,
    messages: list[ChatMessage] | None = None,
    prompt_mode: PromptMode | None = None,
    enable_thinking: bool = True,
    thinking_effort: ThinkingEffort = "medium",
    max_output_tokens: int | None = None,
) -> tuple[int, str]:
    """Call an LLM for single-step retrosynthesis.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.
    model : str, optional
        LiteLLM model identifier, optionally with ``:adv``.
    temperature : float, optional
        Sampling temperature.
    messages : list[ChatMessage], optional
        Explicit chat messages. When provided, default prompt construction is
        skipped.
    prompt_mode : {"standard", "advanced"}, optional
        Prompt-mode override.
    enable_thinking : bool, optional
        Whether provider-supported reasoning controls should be sent.
    thinking_effort : {"low", "medium", "high", "max"}, optional
        Provider reasoning effort when supported.
    max_output_tokens : int, optional
        Output token override. Defaults to the selected prompt configuration.

    Returns
    -------
    tuple[int, str]
        ``(status_code, response_text)``. ``400`` indicates the API call failed
        after retries.

    Examples
    --------
    >>> status, response_text = call_LLM(  # doctest: +SKIP
    ...     "CCO",
    ...     model="openai/gpt-4o-mini",
    ...     max_output_tokens=1024,
    ... )
    >>> isinstance(status, int) and isinstance(response_text, str)  # doctest: +SKIP
    True
    """
    selection = resolve_model_selection(model, prompt_mode=prompt_mode)
    logger.info("Calling LLM", model=selection.completion_model, molecule=molecule)

    addon = build_addon_prompt(molecule)
    default_messages = build_messages(
        molecule=molecule,
        model=model,
        addon=addon,
        messages=messages,
        prompt_mode=prompt_mode,
    )
    _, _, default_max_output_tokens = obtain_prompt(model, prompt_mode=prompt_mode)
    token_limit = max_output_tokens or default_max_output_tokens
    params = build_completion_params(
        model=selection.completion_model,
        messages=default_messages,
        max_completion_tokens=token_limit,
        temperature=temperature,
        enable_thinking=enable_thinking,
        thinking_effort=thinking_effort,
        metadata={"task": "retrosynthesis"},
    )

    last_error = ""
    for attempt in range(1, MAX_API_RETRIES + 1):
        try:
            response = completion(**params)
            content = response.choices[0].message.content
            response_text = coerce_response_text(content)
            logger.debug("Received LLM response", response_length=len(response_text))
            return 200, response_text
        except Exception as exc:
            last_error = str(exc)
            logger.warning(
                "LLM call attempt failed",
                attempt=attempt,
                max_retries=MAX_API_RETRIES,
                model=selection.completion_model,
                error=str(exc),
            )
            if attempt < MAX_API_RETRIES:
                time.sleep(0.5)

    logger.error(
        "All LLM attempts exhausted",
        max_retries=MAX_API_RETRIES,
        model=selection.completion_model,
    )
    return 400, last_error


def parse_cot_response(response_text: str) -> tuple[int, list[str], str]:
    """Parse Claude-style ``<cot>`` plus ``<json>`` output.

    Parameters
    ----------
    response_text : str
        Raw model response containing ``<cot>`` and ``<json>`` tags.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, extracted thinking steps, and JSON payload.

    Examples
    --------
    >>> text = "<cot><thinking>step</thinking></cot><json>{}</json>"
    >>> parse_cot_response(text)
    (200, ['step'], '{}')
    """
    cot_content = extract_tag_content(response_text, "cot")
    json_content = extract_tag_content(response_text, "json")
    if not cot_content or not json_content:
        return 501, [], ""

    thinking_steps = [
        match.group(1)
        for match in re.finditer(
            r"<thinking[^>]*>(.*?)</thinking>",
            cot_content,
            re.DOTALL,
        )
        if match.group(1).strip()
    ]
    if not thinking_steps:
        return 501, [], ""

    return 200, thinking_steps, json_content


def parse_openai_response(response_text: str) -> tuple[int, list[str], str]:
    """Parse OpenAI output into JSON payload content.

    Parameters
    ----------
    response_text : str
        OpenAI response text containing tagged, fenced, or raw JSON.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, empty thinking-step list, and JSON payload.

    Examples
    --------
    >>> parse_openai_response('<json>{"data": []}</json>')
    (200, [], '{"data": []}')
    """
    json_content = extract_json_payload(response_text)
    if not json_content:
        return 502, [], ""
    return 200, [], json_content


def parse_deepseek_response(response_text: str) -> tuple[int, list[str], str]:
    """Parse DeepSeek output into optional thinking and JSON payload content.

    Parameters
    ----------
    response_text : str
        DeepSeek response text containing optional ``<think>`` content and JSON.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, optional thinking steps, and JSON payload.

    Examples
    --------
    >>> parse_deepseek_response('<think>step</think><json>{"data": []}</json>')
    (200, ['step'], '{"data": []}')
    """
    thinking_content = extract_tag_content(response_text, "think")
    json_content = extract_json_payload(response_text)
    if not json_content:
        return 503, [], ""
    thinking_steps = [thinking_content] if thinking_content else []
    return 200, thinking_steps, json_content


def parse_response(response_text: str, model: str) -> tuple[int, list[str], str]:
    """Parse model output into status, visible thinking steps, and JSON content.

    Parameters
    ----------
    response_text : str
        Raw text returned by the LLM.
    model : str
        Model identifier used for the request.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, visible thinking steps, and JSON payload.

    Examples
    --------
    >>> parse_response('<json>{"data": []}</json>', "openai/gpt-4o-mini")
    (200, [], '{"data": []}')
    """
    family = resolve_model_selection(model).family
    try:
        if family == "deepseek":
            return parse_deepseek_response(response_text)
        if family == "openai":
            return parse_openai_response(response_text)
        return parse_cot_response(response_text)
    except Exception as exc:
        logger.error("Error parsing LLM response", family=family, error=str(exc))
        return 505, [], ""


def split_cot_json(response_text: str) -> tuple[int, list[str], str]:
    """Backward-compatible wrapper for parsing Claude-style responses.

    Parameters
    ----------
    response_text : str
        Raw Claude-style response text.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, thinking steps, and JSON payload.

    Examples
    --------
    >>> split_cot_json("<cot><thinking>x</thinking></cot><json>{}</json>")
    (200, ['x'], '{}')
    """
    return parse_cot_response(response_text)


def split_json_openAI(response_text: str) -> tuple[int, str]:
    """Backward-compatible wrapper for parsing OpenAI responses.

    Parameters
    ----------
    response_text : str
        OpenAI response text containing JSON.

    Returns
    -------
    tuple[int, str]
        Status code and JSON payload.

    Examples
    --------
    >>> split_json_openAI('<json>{"data": []}</json>')
    (200, '{"data": []}')
    """
    status, _, json_content = parse_openai_response(response_text)
    return status, json_content


def split_json_deepseek(response_text: str) -> tuple[int, list[str], str]:
    """Backward-compatible wrapper for parsing DeepSeek responses.

    Parameters
    ----------
    response_text : str
        DeepSeek response text containing optional thinking and JSON.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, optional thinking steps, and JSON payload.

    Examples
    --------
    >>> split_json_deepseek('<think>x</think><json>{"data": []}</json>')
    (200, ['x'], '{"data": []}')
    """
    return parse_deepseek_response(response_text)


def split_json_master(response_text: str, model: str) -> tuple[int, list[str], str]:
    """Backward-compatible wrapper for provider-aware response parsing.

    Parameters
    ----------
    response_text : str
        Raw model response text.
    model : str
        Model identifier used to select the parser.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, visible thinking steps, and JSON payload.

    Examples
    --------
    >>> split_json_master('<json>{"data": []}</json>', "openai/gpt-4o-mini")
    (200, [], '{"data": []}')
    """
    return parse_response(response_text, model)


def validate_json_response(
    json_content: str,
) -> tuple[int, list[Pathway], list[str], list[float]]:
    """Parse JSON-like response content into candidate pathways.

    Parameters
    ----------
    json_content : str
        JSON object containing ``data``, ``explanation``, and
        ``confidence_scores``.

    Returns
    -------
    tuple[int, list[Pathway], list[str], list[float]]
        Status code, pathways, explanations, and confidence scores.

    Examples
    --------
    >>> validate_json_response(
    ...     '{"data": [["CCO"]], "explanation": ["demo"], "confidence_scores": [1]}'
    ... )
    (200, [['CCO']], ['demo'], [1.0])
    """
    try:
        result = cast(dict[str, Any], ast.literal_eval(json_content))
        pathways = cast(list[Pathway], result["data"])
        explanations = cast(list[str], result["explanation"])
        confidence = [float(value) for value in result["confidence_scores"]]
        return 200, pathways, explanations, confidence
    except Exception as exc:
        logger.error("Error parsing JSON response", error=str(exc))
        return 504, [], [], []


def validate_split_json(
    json_content: str,
) -> tuple[int, list[Pathway], list[str], list[float]]:
    """Backward-compatible wrapper for response JSON validation.

    Parameters
    ----------
    json_content : str
        JSON object containing pathway data, explanations, and confidence.

    Returns
    -------
    tuple[int, list[Pathway], list[str], list[float]]
        Status code, pathways, explanations, and confidence scores.

    Examples
    --------
    >>> validate_split_json(
    ...     '{"data": [["CCO"]], "explanation": ["demo"], "confidence_scores": [1]}'
    ... )
    (200, [['CCO']], ['demo'], [1.0])
    """
    return validate_json_response(json_content)


def filter_stable_pathways(pathways: list[Pathway]) -> list[Pathway]:
    """Keep pathways whose molecules pass the heuristic stability checker.

    Parameters
    ----------
    pathways : list[Pathway]
        Candidate precursor pathways.

    Returns
    -------
    list[Pathway]
        Pathways whose molecules are assessed as stable enough to keep.

    Examples
    --------
    >>> filter_stable_pathways([["CCO"]])
    [['CCO']]
    """
    stable_pathways: list[Pathway] = []
    for pathway in pathways:
        if all(
            str(check_molecule_stability(smiles).get("assessment", "Unknown"))
            in STABLE_ASSESSMENTS
            for smiles in pathway
        ):
            stable_pathways.append(pathway)
    return stable_pathways


def filter_hallucination_pathways(
    molecule: str,
    pathways: list[Pathway],
) -> list[Pathway]:
    """Keep pathways with low or medium hallucination severity.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.
    pathways : list[Pathway]
        Candidate precursor pathways.

    Returns
    -------
    list[Pathway]
        Pathways whose hallucination severity is acceptable.

    Examples
    --------
    >>> filter_hallucination_pathways("CCO", [["CCO"]])
    [['CCO']]
    """
    filtered_pathways: list[Pathway] = []
    for pathway in pathways:
        combined_smiles = ".".join(pathway)
        report = calculate_hallucination_score(combined_smiles, molecule)
        severity = str(report.get("severity", "critical"))
        if severity in HALLUCINATION_SEVERITIES:
            filtered_pathways.append(pathway)
    return filtered_pathways


def filter_metadata_by_pathways(
    retained_pathways: list[Pathway],
    pathways: list[Pathway],
    explanations: list[str],
    confidence: list[float],
) -> tuple[list[Pathway], list[str], list[float]]:
    """Keep explanations and confidence aligned with filtered pathways.

    Parameters
    ----------
    retained_pathways : list[Pathway]
        Filtered pathways to retain.
    pathways : list[Pathway]
        Original pathway order.
    explanations : list[str]
        Original explanations aligned with ``pathways``.
    confidence : list[float]
        Original confidence scores aligned with ``pathways``.

    Returns
    -------
    tuple[list[Pathway], list[str], list[float]]
        Retained pathways with matching explanations and confidence scores.

    Examples
    --------
    >>> filter_metadata_by_pathways([["B"]], [["A"], ["B"]], ["a", "b"], [0.1, 0.9])
    ([['B']], ['b'], [0.9])
    """
    retained_explanations: list[str] = []
    retained_confidence: list[float] = []
    available = list(zip(pathways, explanations, confidence))

    for retained_pathway in retained_pathways:
        for index, (pathway, explanation, score) in enumerate(available):
            if pathway == retained_pathway:
                retained_explanations.append(explanation)
                retained_confidence.append(score)
                available.pop(index)
                break

    return retained_pathways, retained_explanations, retained_confidence


def llm_pipeline(
    molecule: str,
    model: str = "claude-opus-4-6",
    messages: list[ChatMessage] | None = None,
    stability_check: str | bool = "False",
    hallucination_check: str | bool = "False",
    prompt_mode: PromptMode | None = None,
    enable_thinking: bool = True,
    max_output_tokens: int | None = None,
) -> tuple[list[Pathway], list[str], list[float]]:
    """Run the end-to-end single-step retrosynthesis LLM pipeline.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.
    model : str, optional
        LiteLLM model identifier, optionally with ``:adv``.
    messages : list[ChatMessage], optional
        Explicit chat messages for the call.
    stability_check : str or bool, optional
        Whether to run the heuristic stability filter.
    hallucination_check : str or bool, optional
        Whether to run the heuristic hallucination filter.
    prompt_mode : {"standard", "advanced"}, optional
        Prompt-mode override.
    enable_thinking : bool, optional
        Whether provider-supported reasoning controls should be sent.
    max_output_tokens : int, optional
        Output token override for each LLM call.

    Returns
    -------
    tuple[list[Pathway], list[str], list[float]]
        Valid pathways, explanations, and confidence scores.

    Examples
    --------
    >>> pathways, explanations, confidence = llm_pipeline(  # doctest: +SKIP
    ...     "CCO",
    ...     model="openai/gpt-4o-mini",
    ...     max_output_tokens=2048,
    ... )
    >>> isinstance(pathways, list) and isinstance(explanations, list)  # doctest: +SKIP
    True
    """
    run_stability = is_enabled(stability_check)
    run_hallucination = is_enabled(hallucination_check)
    max_attempts = 15 if (run_stability or run_hallucination) else 6

    for attempt in range(max_attempts):
        temperature = round(attempt * TEMPERATURE_STEP, 2)
        current_model = model
        if resolve_model_selection(model).family == "deepseek" and attempt > 0:
            current_model = FALLBACK_MODEL

        status, response_text = call_LLM(
            molecule=molecule,
            model=current_model,
            messages=messages,
            temperature=temperature,
            prompt_mode=prompt_mode,
            enable_thinking=enable_thinking,
            max_output_tokens=max_output_tokens,
        )
        if status != 200:
            logger.warning("LLM call failed", status_code=status)
            continue

        status, _, json_content = parse_response(response_text, current_model)
        if status != 200:
            logger.warning("Response parsing failed", status_code=status)
            continue

        status, res_molecules, res_explanations, res_confidence = (
            validate_json_response(json_content)
        )
        if status != 200:
            logger.warning("JSON validation failed", status_code=status)
            continue

        output_pathways, output_explanations, output_confidence = validity_check(
            molecule,
            res_molecules,
            res_explanations,
            res_confidence,
        )

        if run_stability and output_pathways:
            previous_pathways = output_pathways
            output_pathways = filter_stable_pathways(previous_pathways)
            output_pathways, output_explanations, output_confidence = (
                filter_metadata_by_pathways(
                    output_pathways,
                    previous_pathways,
                    output_explanations,
                    output_confidence,
                )
            )

        if run_hallucination and output_pathways:
            previous_pathways = output_pathways
            output_pathways = filter_hallucination_pathways(molecule, previous_pathways)
            output_pathways, output_explanations, output_confidence = (
                filter_metadata_by_pathways(
                    output_pathways,
                    previous_pathways,
                    output_explanations,
                    output_confidence,
                )
            )

        if output_pathways:
            logger.info(
                "Pipeline succeeded",
                pathways=len(output_pathways),
                explanations=len(output_explanations),
            )
            return output_pathways, output_explanations, output_confidence

    return [], [], []
