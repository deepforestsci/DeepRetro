"""Helper utilities for DeepRetro LLM provider handling."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Literal, TypedDict

from deepretro.utils.variables import DEEPSEEK_MODELS, OPENAI_MODELS

PromptMode = Literal["standard", "advanced"]
ModelFamily = Literal["deepseek", "openai", "default"]
ProviderName = Literal["anthropic", "openai", "deepseek"]
ThinkingEffort = Literal["low", "medium", "high", "max"]
OutputTokenParam = Literal["max_tokens", "max_completion_tokens"]


class ChatMessage(TypedDict):
    """Chat message sent to LiteLLM."""

    role: str
    content: str


Pathway = list[str]

MIN_REASONING_OUTPUT_TOKENS = 8192
PREFERRED_DEEPSEEK_MODEL = "fireworks_ai/accounts/fireworks/models/deepseek-r1"


@dataclass(frozen=True)
class ModelSelection:
    """Normalized model configuration derived from a model identifier.

    Attributes
    ----------
    raw_model : str
        Original model identifier provided by the caller.
    completion_model : str
        Model identifier passed to LiteLLM after alias normalization.
    prompt_mode : {"standard", "advanced"}
        Prompt variant selected for the call.
    family : {"deepseek", "openai", "default"}
        Prompt and parser family for the model.
    provider : {"anthropic", "openai", "deepseek"}
        Provider inferred from the model name.
    output_token_param : {"max_tokens", "max_completion_tokens"}
        Token-limit keyword expected by the provider.
    supports_reasoning_effort : bool
        Whether the model supports a ``reasoning_effort`` control.
    supports_seed : bool
        Whether deterministic seed should be sent.
    requires_temperature_one : bool
        Whether the model should be called with temperature ``1``.
    """

    raw_model: str
    completion_model: str
    prompt_mode: PromptMode
    family: ModelFamily
    provider: ProviderName
    output_token_param: OutputTokenParam
    supports_reasoning_effort: bool
    supports_seed: bool
    requires_temperature_one: bool


def split_prompt_mode(model: str) -> tuple[str, PromptMode]:
    """Split a model identifier from the optional ``:adv`` prompt suffix.

    Parameters
    ----------
    model : str
        Model identifier, optionally suffixed with ``:adv``.

    Returns
    -------
    tuple[str, PromptMode]
        Base model identifier and prompt mode.

    Examples
    --------
    >>> split_prompt_mode("openai/gpt-4o-mini:adv")
    ('openai/gpt-4o-mini', 'advanced')
    >>> split_prompt_mode("claude-opus-4-6")
    ('claude-opus-4-6', 'standard')
    """
    base_model, separator, suffix = model.rpartition(":")
    if separator and suffix == "adv":
        return base_model, "advanced"
    return model, "standard"


def strip_provider_prefix(model: str) -> str:
    """Return the provider-independent model name when a known prefix exists.

    Parameters
    ----------
    model : str
        Model identifier that may include a LiteLLM provider prefix.

    Returns
    -------
    str
        Model identifier with known ``openai/`` or ``anthropic/`` prefix
        removed.

    Examples
    --------
    >>> strip_provider_prefix("openai/gpt-4o-mini")
    'gpt-4o-mini'
    >>> strip_provider_prefix("fireworks_ai/accounts/fireworks/models/deepseek-r1")
    'fireworks_ai/accounts/fireworks/models/deepseek-r1'
    """
    provider, separator, remainder = model.partition("/")
    if separator and provider in {"openai", "anthropic"}:
        return remainder
    return model


def looks_like_openai_model(model: str) -> bool:
    """Return whether a model identifier appears to be an OpenAI model.

    Parameters
    ----------
    model : str
        Model identifier to inspect.

    Returns
    -------
    bool
        ``True`` for OpenAI-style model identifiers.

    Examples
    --------
    >>> looks_like_openai_model("openai/gpt-4o-mini")
    True
    >>> looks_like_openai_model("claude-opus-4-6")
    False
    """
    base_name = strip_provider_prefix(model).lower()
    return base_name.startswith(("gpt-", "chatgpt-", "o1", "o3", "o4"))


def looks_like_openai_reasoning_model(model: str) -> bool:
    """Return whether a model is an OpenAI reasoning-capable model.

    Parameters
    ----------
    model : str
        Model identifier to inspect.

    Returns
    -------
    bool
        ``True`` for OpenAI reasoning-model families.

    Examples
    --------
    >>> looks_like_openai_reasoning_model("openai/gpt-5")
    True
    >>> looks_like_openai_reasoning_model("openai/gpt-4o-mini")
    False
    """
    base_name = strip_provider_prefix(model).lower()
    return base_name.startswith(("o1", "o3", "o4", "gpt-5"))


def looks_like_anthropic_reasoning_model(model: str) -> bool:
    """Return whether a model is an Anthropic reasoning-capable model.

    Parameters
    ----------
    model : str
        Model identifier to inspect.

    Returns
    -------
    bool
        ``True`` for Anthropic reasoning-model families.

    Examples
    --------
    >>> looks_like_anthropic_reasoning_model("anthropic/claude-sonnet-4-6")
    True
    >>> looks_like_anthropic_reasoning_model("claude-fable-5")
    True
    >>> looks_like_anthropic_reasoning_model("claude-3-5-haiku-20241022")
    False
    """
    base_name = strip_provider_prefix(model).lower()
    # Claude 4.x reasoning families and the Claude 5 family (fable/sonnet-5).
    # These require ``temperature=1`` and accept ``reasoning_effort``.
    return base_name.startswith(
        (
            "claude-opus-4-",
            "claude-sonnet-4-",
            "claude-fable-",
            "claude-sonnet-5",
            "claude-opus-5",
        )
    )


def infer_provider(model: str) -> ProviderName:
    """Infer the provider from a model identifier.

    Parameters
    ----------
    model : str
        Model identifier to classify.

    Returns
    -------
    ProviderName
        Inferred provider name.

    Examples
    --------
    >>> infer_provider("openai/gpt-4o-mini")
    'openai'
    >>> infer_provider("claude-opus-4-6")
    'anthropic'
    """
    lower_model = model.lower()
    if model in DEEPSEEK_MODELS or "deepseek" in lower_model:
        return "deepseek"
    if model in OPENAI_MODELS or lower_model.startswith("openai/"):
        return "openai"
    if looks_like_openai_model(model):
        return "openai"
    if lower_model.startswith("anthropic/") or "claude" in lower_model:
        return "anthropic"
    return "anthropic"


def normalize_completion_model(model: str, provider: ProviderName) -> str:
    """Normalize provider aliases before a LiteLLM completion call.

    Parameters
    ----------
    model : str
        Model identifier supplied by a caller.
    provider : ProviderName
        Provider inferred for the model.

    Returns
    -------
    str
        LiteLLM-compatible model identifier.

    Examples
    --------
    >>> normalize_completion_model("openai/gpt-4o-mini", "openai")
    'openai/gpt-4o-mini'
    >>> normalize_completion_model("fireworks/deepseek-v3p2", "deepseek")
    'fireworks_ai/accounts/fireworks/models/deepseek-r1'
    """
    if provider != "deepseek":
        return model

    lower_model = model.lower()
    if lower_model in {
        "fireworks/deepseek-v3p2",
        "fireworks_ai/deepseek-v3p2",
    }:
        return PREFERRED_DEEPSEEK_MODEL
    if "deepseek" in lower_model and model not in DEEPSEEK_MODELS:
        return PREFERRED_DEEPSEEK_MODEL
    return model


def resolve_model_selection(
    model: str,
    prompt_mode: PromptMode | None = None,
) -> ModelSelection:
    """Resolve provider, prompt, and capability settings for a model.

    Parameters
    ----------
    model : str
        Model identifier, optionally suffixed with ``:adv``.
    prompt_mode : {"standard", "advanced"}, optional
        Explicit prompt-mode override.

    Returns
    -------
    ModelSelection
        Normalized model configuration.

    Examples
    --------
    >>> selection = resolve_model_selection("openai/gpt-4o-mini:adv")
    >>> (selection.provider, selection.prompt_mode)
    ('openai', 'advanced')
    """
    completion_model, suffix_prompt_mode = split_prompt_mode(model)
    provider = infer_provider(completion_model)
    completion_model = normalize_completion_model(completion_model, provider)
    resolved_prompt_mode = prompt_mode or suffix_prompt_mode

    if provider == "deepseek":
        family: ModelFamily = "deepseek"
    elif provider == "openai":
        family = "openai"
    else:
        family = "default"

    is_openai_reasoning = provider == "openai" and looks_like_openai_reasoning_model(
        completion_model
    )
    is_anthropic_reasoning = (
        provider == "anthropic"
        and looks_like_anthropic_reasoning_model(completion_model)
    )
    return ModelSelection(
        raw_model=model,
        completion_model=completion_model,
        prompt_mode=resolved_prompt_mode,
        family=family,
        provider=provider,
        output_token_param=(
            "max_completion_tokens" if provider == "openai" else "max_tokens"
        ),
        supports_reasoning_effort=is_openai_reasoning or is_anthropic_reasoning,
        supports_seed=(provider == "openai") and not is_openai_reasoning,
        requires_temperature_one=is_openai_reasoning or is_anthropic_reasoning,
    )


def resolve_output_token_limit(
    selection: ModelSelection,
    max_output_tokens: int,
    enable_thinking: bool,
) -> int:
    """Return a provider-safe output token budget.

    Parameters
    ----------
    selection : ModelSelection
        Normalized model configuration.
    max_output_tokens : int
        Requested output token limit.
    enable_thinking : bool
        Whether reasoning controls are enabled.

    Returns
    -------
    int
        Token limit adjusted for reasoning-capable models when needed.

    Examples
    --------
    >>> selection = resolve_model_selection("openai/gpt-5")
    >>> resolve_output_token_limit(selection, 32, enable_thinking=True)
    8192
    """
    if enable_thinking and selection.supports_reasoning_effort:
        return max(max_output_tokens, MIN_REASONING_OUTPUT_TOKENS)
    return max_output_tokens


def build_completion_params(
    model: str,
    messages: list[ChatMessage],
    max_completion_tokens: int,
    temperature: float,
    enable_thinking: bool = True,
    thinking_effort: ThinkingEffort = "medium",
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble provider-aware keyword arguments for ``litellm.completion``.

    Parameters
    ----------
    model : str
        Model identifier, optionally with a LiteLLM provider prefix.
    messages : list[ChatMessage]
        Conversation to send to the model.
    max_completion_tokens : int
        Requested maximum output tokens.
    temperature : float
        Sampling temperature for providers that allow it. OpenAI reasoning
        models and Anthropic Claude 4 reasoning-capable models are sent
        ``temperature=1`` when reasoning controls are enabled.
    enable_thinking : bool, optional
        Whether to send reasoning controls for supported models.
    thinking_effort : {"low", "medium", "high", "max"}, optional
        Reasoning effort sent to supported models.
    metadata : dict[str, Any], optional
        Optional LiteLLM metadata for callbacks.

    Returns
    -------
    dict[str, Any]
        Keyword arguments for ``litellm.completion``.

    Examples
    --------
    >>> messages = [{"role": "user", "content": "Reply OK"}]
    >>> params = build_completion_params("openai/gpt-4o-mini", messages, 16, 0.0)
    >>> (params["model"], params["max_completion_tokens"], params["seed"])
    ('openai/gpt-4o-mini', 16, 42)
    >>> params = build_completion_params("anthropic/claude-sonnet-4-6", messages, 16, 0.2)
    >>> (params["max_tokens"], params["temperature"], params["reasoning_effort"])
    (8192, 1, 'medium')
    """
    selection = resolve_model_selection(model)
    output_token_limit = resolve_output_token_limit(
        selection=selection,
        max_output_tokens=max_completion_tokens,
        enable_thinking=enable_thinking,
    )
    params: dict[str, Any] = {
        "model": selection.completion_model,
        "messages": messages,
        selection.output_token_param: output_token_limit,
        "temperature": (
            1 if selection.requires_temperature_one and enable_thinking else temperature
        ),
    }

    if selection.supports_seed:
        params["seed"] = 42
    if metadata is not None:
        params["metadata"] = metadata
    if enable_thinking and selection.supports_reasoning_effort:
        params["reasoning_effort"] = thinking_effort

    return params


def coerce_response_text(content: Any) -> str:
    """Convert LiteLLM response content to a plain string.

    Parameters
    ----------
    content : Any
        Response content returned by a LiteLLM message.

    Returns
    -------
    str
        Plain text response content.

    Examples
    --------
    >>> coerce_response_text("OK")
    'OK'
    >>> coerce_response_text([{"text": "O"}, {"text": "K"}])
    'OK'
    """
    if isinstance(content, str):
        return content
    if content is None:
        return ""
    if isinstance(content, list):
        parts: list[str] = []
        for item in content:
            if isinstance(item, dict):
                text_value = item.get("text")
                if isinstance(text_value, str):
                    parts.append(text_value)
        return "".join(parts)
    return str(content)


def strip_code_fences(text: str) -> str:
    """Remove one surrounding Markdown code fence from a payload.

    Parameters
    ----------
    text : str
        Text that may be wrapped in a Markdown code fence.

    Returns
    -------
    str
        Unwrapped text when a fence is present, otherwise stripped input text.

    Examples
    --------
    >>> strip_code_fences('```json\\n{"data": []}\\n```')
    '{"data": []}'
    >>> strip_code_fences('{"data": []}')
    '{"data": []}'
    """
    stripped = text.strip()
    match = re.fullmatch(r"```(?:json)?\s*(.*?)\s*```", stripped, re.DOTALL)
    return match.group(1).strip() if match else stripped


def extract_tag_content(text: str, tag: str) -> str | None:
    """Extract content between matching XML-like tags.

    Parameters
    ----------
    text : str
        Text containing optional XML-like tags.
    tag : str
        Tag name to extract.

    Returns
    -------
    str | None
        Tag body when found, otherwise ``None``.

    Examples
    --------
    >>> extract_tag_content("<json>{}</json>", "json")
    '{}'
    >>> extract_tag_content("missing", "json") is None
    True
    """
    match = re.search(
        rf"<{re.escape(tag)}\b[^>]*>\s*(.*?)\s*</{re.escape(tag)}>",
        text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def extract_json_payload(response_text: str) -> str | None:
    """Extract tagged, fenced, or raw JSON-like payload from model text.

    Parameters
    ----------
    response_text : str
        Raw model response text.

    Returns
    -------
    str | None
        Extracted JSON-like payload, or ``None`` when no payload is found.

    Examples
    --------
    >>> extract_json_payload('<json>{"data": []}</json>')
    '{"data": []}'
    >>> extract_json_payload('No JSON here') is None
    True
    """
    tagged_json = extract_tag_content(response_text, "json")
    if tagged_json:
        return tagged_json

    stripped = strip_code_fences(response_text)
    if stripped.startswith(("{", "[")) and stripped.endswith(("}", "]")):
        return stripped

    first_object = stripped.find("{")
    last_object = stripped.rfind("}")
    if first_object != -1 and last_object > first_object:
        return stripped[first_object : last_object + 1]

    first_array = stripped.find("[")
    last_array = stripped.rfind("]")
    if first_array != -1 and last_array > first_array:
        return stripped[first_array : last_array + 1]

    return None


def is_enabled(flag: str | bool) -> bool:
    """Normalize string and boolean feature flags.
    
    Note:
    -----
    This function is added for backwards compatibility, ideally we should 
    not have this as type safety is altered with this, we should remove this 
    once we make sure type safety is maintained.

    Parameters
    ----------
    flag : str or bool
        Feature flag value.

    Returns
    -------
    bool
        Normalized boolean value.

    Examples
    --------
    >>> is_enabled("True")
    True
    >>> is_enabled(False)
    False
    """
    if isinstance(flag, bool):
        return flag
    return flag.lower() == "true"
