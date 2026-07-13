"""Single-step tool-calling agent loop for retrosynthesis.

``agentic_single_step`` replaces the non-agentic ``llm_pipeline`` when
``AutoSolver`` runs with ``solve_mode="single_step_agent"``. The model may call
tools (validity / stability / hallucination / ``run_python``) to self-check its
proposed precursors, then emits the same tagged-JSON payload the pipeline
expects. The result is returned in ``llm_pipeline``'s ``(pathways, explanations,
confidence)`` shape and is *still* run through the deterministic safety filters
in :meth:`AutoSolver.run_llm`.

The model call is injectable (``llm_runner``) so the loop is unit-tested with no
network. The default runner calls ``litellm.completion`` with the tool schemas
and normalizes the returned message to an OpenAI-format assistant dict.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from typing import Any, cast

import structlog

from deepretro.agents.tools import build_tool_registry
from deepretro.utils.llm_helpers import ChatMessage, Pathway

logger = structlog.get_logger(__name__)

ModelCall = Callable[[list[dict[str, Any]]], dict[str, Any]]

_TOOL_INSTRUCTION = (
    "\n\nYou may call the provided tools to validate SMILES, check stability, "
    "check for hallucinations, or run Python for any calculation before you "
    "commit to an answer. When you are done, respond with the final answer in "
    "exactly the JSON format described above (do not call a tool in that final "
    "message)."
)


def agentic_single_step(
    molecule: str,
    model: str,
    *,
    tool_backend: str = "structured",
    sandbox: Any | None = None,
    hallucination_checker: Any | None = None,
    max_iterations: int = 12,
    llm_runner: ModelCall | None = None,
    enable_thinking: bool = False,
    max_output_tokens: int | None = None,
) -> tuple[list[Pathway], list[str], list[float]]:
    """Propose precursors for one molecule via a tool-calling agent.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.
    model : str
        LiteLLM model identifier.
    tool_backend : {"structured", "sandbox"}, optional
        Which tools to expose. ``sandbox`` adds ``run_python``.
    sandbox : Sandbox or None, optional
        Sandbox for ``run_python`` (a default is created if needed).
    hallucination_checker : callable or None, optional
        Resolved checker exposed as the ``check_hallucination`` tool.
    max_iterations : int, optional
        Maximum model turns before giving up.
    llm_runner : callable, optional
        Injectable model call ``(messages) -> assistant_message_dict``. When
        ``None``, ``litellm.completion`` is used.
    enable_thinking : bool, optional
        Reasoning controls. Defaults to ``False``: multi-turn tool use with
        extended thinking is rejected by Anthropic unless the signed thinking
        block is preserved.
    max_output_tokens : int, optional
        Output-token override for the model call.

    Returns
    -------
    tuple[list[Pathway], list[str], list[float]]
        Pathways, explanations, and confidence scores. Empty lists when the
        agent produces no parseable final answer within ``max_iterations``.

    Examples
    --------
    >>> final = {
    ...     "role": "assistant",
    ...     "content": '<json>{"data": [["CCO"]], "explanation": ["reduce"], '
    ...     '"confidence_scores": [0.8]}</json>',
    ... }
    >>> agentic_single_step("CC=O", "openai/gpt-4o-mini", llm_runner=lambda m: final)
    ([['CCO']], ['reduce'], [0.8])
    """
    registry = build_tool_registry(hallucination_checker, sandbox, tool_backend)
    messages = _build_initial_messages(molecule, model)
    call_model = llm_runner or _make_default_model_call(
        model, registry.schemas, enable_thinking, max_output_tokens
    )

    for _iteration in range(max_iterations):
        assistant = call_model(messages)
        messages.append(assistant)

        tool_calls = assistant.get("tool_calls")
        if not tool_calls:
            return _parse_final_answer(assistant.get("content") or "", model)

        for tool_call in tool_calls:
            function = tool_call.get("function", {})
            name = function.get("name", "")
            arguments = _parse_arguments(function.get("arguments"))
            result = registry.execute(name, arguments)
            logger.info("Agent tool call", tool=name, arguments=arguments, result=result)
            messages.append(
                {
                    "role": "tool",
                    "tool_call_id": tool_call.get("id", ""),
                    "content": json.dumps(result),
                }
            )

    logger.warning(
        "Agent reached max_iterations without a final answer",
        molecule=molecule,
        max_iterations=max_iterations,
    )
    return [], [], []


def agentic_orchestrator(
    molecule: str,
    model: str,
    **kwargs: Any,
) -> tuple[list[Pathway], list[str], list[float]]:
    """Top-level tool-driven search over the whole tree (not implemented).

    The interface and flag wiring exist; ``single_step_agent`` is the built-out
    agent mode. This scaffold raises so the reserved mode fails loudly.

    Raises
    ------
    NotImplementedError
        Always.

    Examples
    --------
    >>> agentic_orchestrator("CCO", "openai/gpt-4o-mini")
    Traceback (most recent call last):
    NotImplementedError: agentic_orchestrator is not implemented yet
    """
    raise NotImplementedError("agentic_orchestrator is not implemented yet")


def _build_initial_messages(molecule: str, model: str) -> list[dict[str, Any]]:
    """Build the [system, user] messages with tool-usage guidance appended."""
    from deepretro.utils.llm import build_messages

    messages = [dict(message) for message in build_messages(molecule, model)]
    if messages and messages[0].get("role") == "system":
        messages[0]["content"] = f"{messages[0].get('content', '')}{_TOOL_INSTRUCTION}"
    return messages


def _parse_arguments(arguments: Any) -> dict[str, Any]:
    """Parse tool-call arguments (a JSON string or already a dict)."""
    if isinstance(arguments, dict):
        return arguments
    if not arguments:
        return {}
    try:
        parsed = json.loads(arguments)
    except (json.JSONDecodeError, TypeError):
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _parse_final_answer(
    content: str,
    model: str,
) -> tuple[list[Pathway], list[str], list[float]]:
    """Parse a final assistant message into pathways/explanations/confidence."""
    from deepretro.utils.llm import parse_response, validate_split_json

    status, _thinking, json_content = parse_response(content, model)
    if status != 200 or not json_content:
        return [], [], []
    status, pathways, explanations, confidence = validate_split_json(json_content)
    if status != 200:
        return [], [], []
    return pathways, explanations, confidence


def _make_default_model_call(
    model: str,
    tools: list[dict[str, Any]],
    enable_thinking: bool,
    max_output_tokens: int | None,
) -> ModelCall:
    """Build the default ``litellm.completion``-backed model call."""

    def _call(messages: list[dict[str, Any]]) -> dict[str, Any]:
        from litellm import completion

        from deepretro.utils.llm_helpers import build_completion_params

        params = build_completion_params(
            model=model,
            # The conversation is an OpenAI-format superset of ChatMessage
            # (assistant tool-call turns, tool results); litellm accepts it.
            messages=cast("list[ChatMessage]", messages),
            max_completion_tokens=max_output_tokens or 8192,
            temperature=0.0,
            enable_thinking=enable_thinking,
            metadata={"task": "retrosynthesis_agent"},
        )
        params["tools"] = tools
        response = completion(**params)
        message = response.choices[0].message
        return message.model_dump()

    return _call
