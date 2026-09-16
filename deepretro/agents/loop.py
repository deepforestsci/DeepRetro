"""Single-step tool-calling agent loop for retrosynthesis.

``agentic_single_step`` replaces the non-agentic ``llm_pipeline`` when
``AutoSolver`` runs with ``solve_mode="single_step_agent"``. The model may call
tools (validity / stability / hallucination / ``run_python``) to self-check its
proposed precursors, then emits the same tagged-JSON payload the pipeline
expects. The result is returned in ``llm_pipeline``'s ``(pathways, explanations,
confidence)`` shape and is *still* run through the deterministic safety filters
in :meth:`AutoSolver.run_llm`.

The model call is injectable (``llm_runner``) so the loop is unit-tested with no
network. Every returned assistant message crosses a lossless history boundary
before the loop inspects tool calls or visible content.
"""

from __future__ import annotations

import json
import math
import time
from collections.abc import Callable
from typing import Any, cast

import structlog

from deepretro.agents.message_history import append_assistant_message
from deepretro.agents.tools import build_tool_registry
from deepretro.utils.llm_helpers import ChatMessage, Pathway
from deepretro.utils.llm_trace import elapsed_ms, langfuse_metadata, record_llm_call

logger = structlog.get_logger(__name__)

ModelCall = Callable[[list[dict[str, Any]]], Any]

DEFAULT_MIN_ITERATIONS = 5
DEFAULT_MAX_ITERATIONS = 15
DEFAULT_ITERATION_DECAY = 0.75


def validate_iteration_budget(
    min_iterations: int, max_iterations: int, decay: float
) -> None:
    """Validate the parameters of :func:`iteration_budget`.

    Parameters
    ----------
    min_iterations : int
        Lower cap on the budget; must be at least 1.
    max_iterations : int
        Upper cap on the budget; must be at least ``min_iterations``.
    decay : float
        Per-depth multiplier; must satisfy ``0 < decay <= 1``.

    Raises
    ------
    ValueError
        If any parameter is out of range.

    Examples
    --------
    >>> validate_iteration_budget(5, 15, 0.75)
    >>> validate_iteration_budget(6, 5, 0.75)
    Traceback (most recent call last):
        ...
    ValueError: max_iterations must be >= min_iterations
    """
    if min_iterations < 1:
        raise ValueError("min_iterations must be at least 1")
    if max_iterations < min_iterations:
        raise ValueError("max_iterations must be >= min_iterations")
    if not 0.0 < decay <= 1.0:
        raise ValueError("decay must be in the interval (0, 1]")


def iteration_budget(
    smiles: str,
    depth: int,
    *,
    min_iterations: int = DEFAULT_MIN_ITERATIONS,
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    decay: float = DEFAULT_ITERATION_DECAY,
) -> int:
    """Size the agent's turn budget for one node from carbon count and depth.

    The budget starts at the number of carbon atoms in ``smiles`` and shrinks
    geometrically with recursion depth (``carbons * decay ** depth``), so a
    large target gets many turns at the root while small, deep intermediates
    get few. The value is rounded half-up and clamped to
    ``[min_iterations, max_iterations]``. Molecules RDKit cannot parse, or
    that contain no carbon, receive ``min_iterations``.

    Parameters
    ----------
    smiles : str
        Molecule at the current node.
    depth : int
        Recursion depth of the node (``0`` for the target).
    min_iterations : int, optional
        Lower cap on the budget. Defaults to 5.
    max_iterations : int, optional
        Upper cap on the budget. Defaults to 15.
    decay : float, optional
        Multiplier applied once per depth level. Defaults to 0.75.

    Returns
    -------
    int
        Maximum number of model turns for this node.

    Raises
    ------
    ValueError
        If ``depth`` is negative or the caps/decay are out of range.

    Examples
    --------
    >>> iteration_budget("C" * 20, depth=0)
    15
    >>> iteration_budget("C" * 20, depth=2)
    11
    >>> iteration_budget("CCO", depth=0)
    5
    """
    from rdkit import Chem

    validate_iteration_budget(min_iterations, max_iterations, decay)
    if depth < 0:
        raise ValueError("depth must be non-negative")

    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return min_iterations
    carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    raw = math.floor(carbons * decay**depth + 0.5)
    return max(min_iterations, min(max_iterations, raw))


_TOOL_INSTRUCTION = (
    "\n\nYou may call the provided tools to validate SMILES, check stability, "
    "check for hallucinations, or run Python for any calculation before you "
    "commit to an answer. When existing protecting groups distract from the "
    "core transformation, call handle_protection to mask supported motifs. "
    "Treat the returned mapped dummy atoms as unchanged protected substituents, "
    "not disconnection targets. Motif matches alone do not establish a protecting "
    "role; consider compatibility with the proposed conditions. After reasoning "
    "about the exposed core, call handle_deprotection with the returned mask_id "
    "and mode='restore' for each masked precursor to restore the full structure. "
    "To explore actual chemical deprotection, use handle_deprotection with "
    "mode='propose' and full protected SMILES, without a mask_id. Those candidates "
    "describe forward protected-substrate to deprotected-product transformations, "
    "not retro precursors of the protected input. For a retro deprotection step, "
    "test your proposed protected precursor and compare its deprotected product "
    "with the target. Assess conditions, selectivity, "
    "and compatibility; the tool does not validate these. Validate restored SMILES "
    "and never include dummy atoms in final pathways. When you are done, "
    "respond with the final answer in "
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
    max_iterations: int = 6,
    llm_runner: ModelCall | None = None,
    enable_thinking: bool = True,
    max_output_tokens: int | None = None,
    event_sink: list[dict[str, Any]] | None = None,
    az_tools: bool = False,
    az_model: str = "Pistachio_100+",
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
        Injectable model call ``(messages) -> assistant_message``. The result
        may be a serialized mapping or an object exposing ``model_dump()``.
        When ``None``, ``litellm.completion`` is used.
    enable_thinking : bool, optional
        Whether provider-supported reasoning controls should be enabled.
        Defaults to ``True``. Provider-returned reasoning artifacts are treated
        as opaque protocol data and retained in memory only for this agent run
        so signed multi-turn tool conversations can be replayed losslessly.
        Pass ``False`` to disable provider reasoning explicitly.
    max_output_tokens : int, optional
        Output-token override for the model call.
    event_sink : list of dict or None, optional
        If provided, agent turns that yield no usable pathway append an event
        describing what happened: a ``refusal`` (the model declined), a
        ``no_parseable_answer`` (a final message that did not parse), or a
        ``no_final_answer`` (``max_iterations`` reached without a final
        message). Lets callers surface model refusals in their output.

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
    registry = build_tool_registry(
        hallucination_checker,
        sandbox,
        tool_backend,
        az_tools=az_tools,
        az_model=az_model,
    )
    messages = _build_initial_messages(molecule, model)
    call_model = llm_runner or _make_default_model_call(
        model, registry.schemas, enable_thinking, max_output_tokens
    )

    for _iteration in range(max_iterations):
        assistant = append_assistant_message(messages, call_model(messages))

        tool_calls = assistant.get("tool_calls")
        if not tool_calls:
            content = assistant.get("content") or ""
            result = _parse_final_answer(content, model)
            if _contains_dummy_atoms(result[0]):
                messages.append(
                    {
                        "role": "user",
                        "content": (
                            "Your answer contains unresolved dummy atoms. Restore each "
                            "masked precursor with handle_deprotection and its mask_id, "
                            "then return full molecular SMILES in the required JSON format."
                        ),
                    }
                )
                continue
            if event_sink is not None and not result[0]:
                event_sink.append(_classify_agent_event(molecule, content))
            return result

        for tool_call in tool_calls:
            function = tool_call.get("function", {})
            name = function.get("name", "")
            arguments = _parse_arguments(function.get("arguments"))
            result = registry.execute(name, arguments)
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
    if event_sink is not None:
        event_sink.append(
            {
                "molecule": molecule,
                "kind": "no_final_answer",
                "detail": f"reached max_iterations={max_iterations}",
            }
        )
    return [], [], []


def _contains_dummy_atoms(pathways: list[Pathway]) -> bool:
    """Check parsed pathways for graph placeholders, which are not molecules."""
    from rdkit import Chem

    for pathway in pathways:
        for smiles in pathway:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None and any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                return True
    return False


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
    """Build the system and user messages with tool guidance appended.

    Examples
    --------
    >>> messages = _build_initial_messages("CCO", "openai/gpt-4o-mini")
    >>> [message["role"] for message in messages]
    ['system', 'user']
    >>> "You may call" in messages[0]["content"]
    True
    """
    from deepretro.utils.llm import build_messages

    messages = [dict(message) for message in build_messages(molecule, model)]
    if messages and messages[0].get("role") == "system":
        messages[0]["content"] = f"{messages[0].get('content', '')}{_TOOL_INSTRUCTION}"
    return messages


_REFUSAL_MARKERS = (
    "i can't",
    "i cannot",
    "i can not",
    "i'm unable",
    "i am unable",
    "i won't",
    "i will not",
    "cannot assist",
    "can't help",
    "cannot help",
    "not able to",
    "i must decline",
    "i refuse",
    "against my",
    "unable to provide",
    "cannot provide",
    "can't provide",
    "not comfortable",
)


def _classify_agent_event(molecule: str, content: str) -> dict[str, Any]:
    """Classify an empty agent final message as a refusal vs unparseable answer.

    Examples
    --------
    >>> _classify_agent_event("CCO", "I can't help with that.")["kind"]
    'refusal'
    >>> _classify_agent_event("CCO", "here is my answer")["kind"]
    'no_parseable_answer'
    """
    lowered = content.lower()
    kind = (
        "refusal"
        if any(m in lowered for m in _REFUSAL_MARKERS)
        else "no_parseable_answer"
    )
    return {"molecule": molecule, "kind": kind, "content_excerpt": content[:2000]}


def _parse_arguments(arguments: Any) -> dict[str, Any]:
    """Parse tool-call arguments from a JSON string or existing dictionary.

    Examples
    --------
    >>> _parse_arguments('{"smiles": "CCO"}')
    {'smiles': 'CCO'}
    >>> _parse_arguments("not JSON")
    {}
    """
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
    """Parse a final assistant message into pipeline result lists.

    Examples
    --------
    >>> content = (
    ...     '<json>{"data": [["CCO"]], "explanation": ["reduce"], '
    ...     '"confidence_scores": [0.8]}</json>'
    ... )
    >>> _parse_final_answer(content, "openai/gpt-4o-mini")
    ([['CCO']], ['reduce'], [0.8])
    """
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
    """Build the default ``litellm.completion``-backed model call.

    Examples
    --------
    >>> call_model = _make_default_model_call(
    ...     "openai/gpt-4o-mini", [], True, 1024
    ... )
    >>> callable(call_model)
    True
    """

    # One model call per agent loop turn, so the call count is the iteration.
    iteration = 0

    def _call(messages: list[dict[str, Any]]) -> object:
        """Send one conversation turn to LiteLLM and return its raw message.

        Examples
        --------
        The closure is obtained through ``_make_default_model_call``:

        >>> call_model = _make_default_model_call(
        ...     "openai/gpt-4o-mini", [], True, 1024
        ... )
        >>> callable(call_model)
        True
        """
        nonlocal iteration

        from litellm import completion

        from deepretro.utils.llm_helpers import build_completion_params

        iteration += 1
        params = build_completion_params(
            model=model,
            # The conversation is an OpenAI-format superset of ChatMessage
            # (assistant tool-call turns, tool results); litellm accepts it.
            messages=cast("list[ChatMessage]", messages),
            max_completion_tokens=max_output_tokens or 8192,
            temperature=0.0,
            enable_thinking=enable_thinking,
            metadata=langfuse_metadata(
                {"task": "retrosynthesis_agent"}, stage="retrosynthesis_agent"
            ),
        )
        params["tools"] = tools
        started = time.perf_counter()
        try:
            response = completion(**params)
        except Exception as exc:
            record_llm_call(
                stage="retrosynthesis_agent",
                model=model,
                messages=messages,
                response=None,
                error=str(exc),
                latency_ms=elapsed_ms(started),
                iteration=iteration,
            )
            raise
        latency_ms = elapsed_ms(started)
        message = response.choices[0].message
        record_llm_call(
            stage="retrosynthesis_agent",
            model=model,
            messages=messages,
            response=response,
            latency_ms=latency_ms,
            tool_calls=getattr(message, "tool_calls", None),
            iteration=iteration,
        )
        return message

    return _call
