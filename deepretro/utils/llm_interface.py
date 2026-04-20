"""Provider-specific LLM interfaces for DeepRetro retrosynthesis calls."""

from __future__ import annotations

import re
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass

import structlog
from litellm import completion

from deepretro.utils.llm_helpers import (
    ChatMessage,
    ModelSelection,
    PromptMode,
    ThinkingEffort,
    build_completion_params,
    coerce_response_text,
    extract_json_payload,
    extract_tag_content,
    resolve_model_selection,
)
from deepretro.utils.utils_molecule import detect_seven_member_rings
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

logger = structlog.get_logger(__name__)

MAX_API_RETRIES = 2
PROMPT_CONFIG = {
    ("deepseek", "standard"): (SYS_PROMPT_DEEPSEEK, USER_PROMPT_DEEPSEEK, 16384),
    ("deepseek", "advanced"): (SYS_PROMPT_V4, USER_PROMPT_DEEPSEEK_V4, 16384),
    ("openai", "standard"): (SYS_PROMPT_OPENAI, USER_PROMPT_OPENAI, 16384),
    ("openai", "advanced"): (SYS_PROMPT_OPENAI, USER_PROMPT_OPENAI, 16384),
    ("default", "standard"): (SYS_PROMPT, USER_PROMPT, 16384),
    ("default", "advanced"): (SYS_PROMPT_V4, USER_PROMPT_V4, 16384),
}


@dataclass(frozen=True)
class LLMRequest:
    """Input state for one LLM completion request.

    Attributes
    ----------
    molecule : str
        Target molecule SMILES.
    model : str
        LiteLLM model identifier, optionally suffixed with ``:adv``.
    temperature : float, optional
        Sampling temperature.
    messages : list[ChatMessage], optional
        Explicit chat messages. When provided, prompt construction is skipped.
    prompt_mode : {"standard", "advanced"}, optional
        Prompt-mode override.
    enable_thinking : bool, optional
        Whether provider-supported reasoning controls should be sent.
    thinking_effort : {"low", "medium", "high", "max"}, optional
        Provider reasoning effort when supported.
    max_output_tokens : int, optional
        Output token override.

    Examples
    --------
    >>> request = LLMRequest("CCO", "openai/gpt-4o-mini", max_output_tokens=16)
    >>> (request.molecule, request.temperature, request.enable_thinking)
    ('CCO', 0.0, True)
    """

    molecule: str
    model: str
    temperature: float = 0.0
    messages: list[ChatMessage] | None = None
    prompt_mode: PromptMode | None = None
    enable_thinking: bool = True
    thinking_effort: ThinkingEffort = "medium"
    max_output_tokens: int | None = None


@dataclass(frozen=True)
class LLMResponse:
    """Output state from one LLM completion request.

    Attributes
    ----------
    status_code : int
        Status code for the call. ``200`` means success.
    text : str
        Response text or final error message.

    Examples
    --------
    >>> response = LLMResponse(status_code=200, text="OK")
    >>> (response.status_code, response.text)
    (200, 'OK')
    """

    status_code: int
    text: str


def build_addon_prompt(molecule: str) -> str:
    """Build additional prompt context for a target molecule.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.

    Returns
    -------
    str
        Additional prompt text for known chemistry edge cases.

    Examples
    --------
    >>> build_addon_prompt("CCO")
    ''
    """
    if detect_seven_member_rings(molecule):
        logger.debug("Detected seven-membered ring", molecule=molecule)
        return ADDON_PROMPT_7_MEMBER
    return ""


class LLMInterface(ABC):
    """Abstract provider interface for LLM completion calls.

    Parameters
    ----------
    model : str
        LiteLLM model identifier.
    prompt_mode : {"standard", "advanced"}, optional
        Prompt-mode override.

    Examples
    --------
    >>> interface = create_llm_interface("openai/gpt-4o-mini")
    >>> isinstance(interface, LLMInterface)
    True
    """

    def __init__(self, model: str, prompt_mode: PromptMode | None = None) -> None:
        """Initialize provider selection for a model."""
        self.selection = resolve_model_selection(model, prompt_mode=prompt_mode)

    def obtain_prompt(self) -> tuple[str, str, int]:
        """Return the prompt pair and default output-token limit.

        Returns
        -------
        tuple[str, str, int]
            System prompt, user prompt, and default output-token limit.

        Examples
        --------
        >>> create_llm_interface("openai/gpt-4o-mini").obtain_prompt()[2]
        16384
        """
        return PROMPT_CONFIG[(self.selection.family, self.selection.prompt_mode)]

    def build_messages(self, request: LLMRequest) -> list[ChatMessage]:
        """Return caller messages or build the default prompt pair.

        Parameters
        ----------
        request : LLMRequest
            Request state for the call.

        Returns
        -------
        list[ChatMessage]
            Chat messages ready for LiteLLM.

        Examples
        --------
        >>> request = LLMRequest("CCO", "openai/gpt-4o-mini")
        >>> messages = create_llm_interface(request.model).build_messages(request)
        >>> [message["role"] for message in messages]
        ['system', 'user']
        """
        if request.messages is not None:
            return request.messages

        sys_prompt, user_prompt, _ = self.obtain_prompt()
        addon = build_addon_prompt(request.molecule)
        user_content = user_prompt.replace("{target_smiles}", request.molecule)
        if addon:
            user_content = f"{user_content}\n\n{addon}"

        return [
            {"role": "system", "content": sys_prompt + addon},
            {"role": "user", "content": user_content},
        ]

    def build_completion_params(
        self,
        request: LLMRequest,
        messages: list[ChatMessage],
    ) -> dict[str, object]:
        """Build provider-aware LiteLLM completion parameters.

        Parameters
        ----------
        request : LLMRequest
            Request state for the call.
        messages : list[ChatMessage]
            Chat messages to send.

        Returns
        -------
        dict[str, object]
            Keyword arguments for ``litellm.completion``.

        Examples
        --------
        >>> request = LLMRequest("CCO", "openai/gpt-4o-mini", max_output_tokens=16)
        >>> params = create_llm_interface(request.model).build_completion_params(
        ...     request,
        ...     [{"role": "user", "content": "Reply OK"}],
        ... )
        >>> params["max_completion_tokens"]
        16
        """
        _, _, default_max_output_tokens = self.obtain_prompt()
        token_limit = request.max_output_tokens or default_max_output_tokens
        return build_completion_params(
            model=self.selection.completion_model,
            messages=messages,
            max_completion_tokens=token_limit,
            temperature=request.temperature,
            enable_thinking=request.enable_thinking,
            thinking_effort=request.thinking_effort,
            metadata={"task": "retrosynthesis"},
        )

    @abstractmethod
    def parse_response(self, response_text: str) -> tuple[int, list[str], str]:
        """Parse provider response text into thinking steps and JSON content.

        Parameters
        ----------
        response_text : str
            Raw model response text.

        Returns
        -------
        tuple[int, list[str], str]
            Status code, visible thinking steps, and JSON payload.

        Examples
        --------
        >>> create_llm_interface("openai/gpt-4o-mini").parse_response(
        ...     '<json>{"data": []}</json>'
        ... )
        (200, [], '{"data": []}')
        """

    def call(self, request: LLMRequest) -> LLMResponse:
        """Call the provider through LiteLLM.

        Parameters
        ----------
        request : LLMRequest
            Request state for the call.

        Returns
        -------
        LLMResponse
            Status code and response text.

        Examples
        --------
        >>> response = create_llm_interface("openai/gpt-4o-mini").call(  # doctest: +SKIP
        ...     LLMRequest("CCO", "openai/gpt-4o-mini")
        ... )
        >>> isinstance(response.text, str)  # doctest: +SKIP
        True
        """
        logger.info(
            "Calling LLM",
            model=self.selection.completion_model,
            molecule=request.molecule,
        )
        messages = self.build_messages(request)
        params = self.build_completion_params(request, messages)

        last_error = ""
        for attempt in range(1, MAX_API_RETRIES + 1):
            try:
                response = completion(**params)
                content = response.choices[0].message.content
                response_text = coerce_response_text(content)
                logger.debug(
                    "Received LLM response", response_length=len(response_text)
                )
                return LLMResponse(status_code=200, text=response_text)
            except Exception as exc:
                last_error = str(exc)
                logger.warning(
                    "LLM call attempt failed",
                    attempt=attempt,
                    max_retries=MAX_API_RETRIES,
                    model=self.selection.completion_model,
                    error=str(exc),
                )
                if attempt < MAX_API_RETRIES:
                    time.sleep(0.5)

        logger.error(
            "All LLM attempts exhausted",
            max_retries=MAX_API_RETRIES,
            model=self.selection.completion_model,
        )
        return LLMResponse(status_code=400, text=last_error)


class AnthropicLLM(LLMInterface):
    """LLM interface for Anthropic/Claude-compatible responses.

    Examples
    --------
    >>> AnthropicLLM("claude-opus-4-6").parse_response(
    ...     "<cot><thinking>x</thinking></cot><json>{}</json>"
    ... )
    (200, ['x'], '{}')
    """

    def parse_response(self, response_text: str) -> tuple[int, list[str], str]:
        """Parse Claude-style ``<cot>`` output with a JSON payload.

        Parameters
        ----------
        response_text : str
            Raw Claude-style response text. The response must contain ``<cot>``
            with one or more ``<thinking>`` entries. JSON may be in ``<json>``
            tags or in a fenced/raw JSON payload after ``</cot>``.

        Returns
        -------
        tuple[int, list[str], str]
            Status code, thinking steps, and JSON payload.

        Examples
        --------
        >>> AnthropicLLM("claude-opus-4-6").parse_response(
        ...     "<cot><thinking>x</thinking></cot><json>{}</json>"
        ... )
        (200, ['x'], '{}')
        >>> AnthropicLLM("claude-opus-4-6").parse_response(
        ...     '<cot><thinking>x</thinking></cot>```json\\n{}\\n```'
        ... )
        (200, ['x'], '{}')
        """
        return parse_cot_response(response_text)


class OpenAILLM(LLMInterface):
    """LLM interface for OpenAI responses.

    Examples
    --------
    >>> OpenAILLM("openai/gpt-4o-mini").parse_response('<json>{"data": []}</json>')
    (200, [], '{"data": []}')
    """

    def parse_response(self, response_text: str) -> tuple[int, list[str], str]:
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
        >>> OpenAILLM("openai/gpt-4o-mini").parse_response('<json>{"data": []}</json>')
        (200, [], '{"data": []}')
        """
        json_content = extract_json_payload(response_text)
        if not json_content:
            return 502, [], ""
        return 200, [], json_content


class DeepSeekLLM(LLMInterface):
    """LLM interface for DeepSeek responses.

    Examples
    --------
    >>> DeepSeekLLM("deepseek-ai/DeepSeek-R1").parse_response(
    ...     '<think>x</think><json>{"data": []}</json>'
    ... )
    (200, ['x'], '{"data": []}')
    """

    def parse_response(self, response_text: str) -> tuple[int, list[str], str]:
        """Parse DeepSeek output into optional thinking and JSON payload content.

        Parameters
        ----------
        response_text : str
            DeepSeek response text containing optional ``<think>`` content and
            tagged, fenced, or raw JSON.

        Returns
        -------
        tuple[int, list[str], str]
            Status code, optional thinking steps, and JSON payload.

        Examples
        --------
        >>> DeepSeekLLM("deepseek-ai/DeepSeek-R1").parse_response(
        ...     '<think>x</think><json>{"data": []}</json>'
        ... )
        (200, ['x'], '{"data": []}')
        """
        thinking_content = extract_tag_content(response_text, "think")
        json_content = extract_json_payload(response_text)
        if not json_content:
            return 503, [], ""
        thinking_steps = [thinking_content] if thinking_content else []
        return 200, thinking_steps, json_content





def parse_cot_response(response_text: str) -> tuple[int, list[str], str]:
    """Parse Claude-style ``<cot>`` output with a JSON payload.

    Parameters
    ----------
    response_text : str
        Raw model response containing ``<cot>`` thinking content and either a
        ``<json>`` payload or a fenced/raw JSON payload after ``</cot>``.

    Returns
    -------
    tuple[int, list[str], str]
        Status code, extracted thinking steps, and JSON payload.

    Examples
    --------
    >>> parse_cot_response("<cot><thinking>x</thinking></cot><json>{}</json>")
    (200, ['x'], '{}')
    >>> parse_cot_response('<cot><thinking>x</thinking></cot>```json\\n{}\\n```')
    (200, ['x'], '{}')
    """
    cot_content = extract_tag_content(response_text, "cot")
    json_content = extract_tag_content(response_text, "json")
    if not json_content:
        response_after_cot = response_text.split("</cot>", maxsplit=1)[-1]
        json_content = extract_json_payload(response_after_cot)
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


def create_llm_interface(
    model: str,
    prompt_mode: PromptMode | None = None,
) -> LLMInterface:
    """Create a provider-specific LLM interface.

    Parameters
    ----------
    model : str
        LiteLLM model identifier.
    prompt_mode : {"standard", "advanced"}, optional
        Prompt-mode override.

    Returns
    -------
    LLMInterface
        Provider-specific interface instance.

    Examples
    --------
    >>> type(create_llm_interface("openai/gpt-4o-mini")).__name__
    'OpenAILLM'
    """
    selection: ModelSelection = resolve_model_selection(model, prompt_mode=prompt_mode)
    if selection.provider == "openai":
        return OpenAILLM(model, prompt_mode=prompt_mode)
    if selection.provider == "deepseek":
        return DeepSeekLLM(model, prompt_mode=prompt_mode)
    if selection.provider == "anthropic":
        return AnthropicLLM(model, prompt_mode=prompt_mode)
    return AnthropicLLM(model, prompt_mode=prompt_mode)
