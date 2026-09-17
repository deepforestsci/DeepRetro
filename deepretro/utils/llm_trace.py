"""Per-molecule LLM trace context for Langfuse grouping and local call logs.

Every LLM call DeepRetro makes while solving one target molecule belongs to the
same logical unit of work, but LiteLLM's Langfuse callback only sees one call at
a time. This module carries that missing grouping in a
:class:`contextvars.ContextVar`:

* :func:`molecule_trace` opens a trace for one target molecule and gives it a
  stable ``session_id``. :func:`langfuse_metadata` injects that id (plus the
  trace / generation names LiteLLM's Langfuse integration understands) into the
  ``metadata`` of every completion call, so Langfuse groups the whole
  retrosynthesis under one session.
* :func:`record_llm_call` mirrors each call to a local ``llm_calls.jsonl`` file
  inside the molecule's output directory, so a run can be inspected offline and
  without a Langfuse account.
* :func:`record_tool_results` logs every tool the agent ran (name, arguments,
  output) both to that file and, when Langfuse credentials are configured, as
  Langfuse *events* inside the same trace, so tool outputs render online next
  to the model turns that requested them.

Both are best-effort observability: neither ever raises into the solver.

Examples
--------
>>> with molecule_trace("CCO", session_id="demo") as trace:
...     langfuse_metadata({"task": "retrosynthesis"}, stage="retrosynthesis")[
...         "session_id"
...     ]
'demo'
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
import uuid
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import structlog

logger = structlog.get_logger(__name__)

#: Name of the per-molecule JSONL log written next to ``pathway_<i>.json``.
LOG_FILENAME = "llm_calls.jsonl"

#: Langfuse trace name shared by every generation of one retrosynthesis run.
TRACE_NAME = "autosolve"

_trace_var: ContextVar[MoleculeTrace | None] = ContextVar(
    "deepretro_molecule_trace", default=None
)
_depth_var: ContextVar[int] = ContextVar("deepretro_node_depth", default=0)
_node_molecule_var: ContextVar[str | None] = ContextVar(
    "deepretro_node_molecule", default=None
)

# A broken log destination must warn once, not once per LLM call.
_write_failure_logged = False

#: Message keys that hold opaque provider replay data (Anthropic thinking
#: signatures). Required on the wire, useless and large in a log.
_OPAQUE_MESSAGE_KEYS = ("thinking_blocks", "provider_specific_fields")

_client_cache: Any | None = None
_client_checked = False


@dataclass(frozen=True)
class MoleculeTrace:
    """Grouping context for every LLM call made for one target molecule.

    Parameters
    ----------
    molecule : str
        Target (root) molecule SMILES the current run is solving.
    session_id : str
        Langfuse session id shared by every generation of this run.
    log_path : Path or None
        Destination of the local JSONL call log, or ``None`` when local
        logging is disabled.
    trace_id : str
        Langfuse trace id shared by every generation and tool event of this
        run. Passed to LiteLLM as ``metadata["trace_id"]`` and used directly
        by the Langfuse client for tool events. Generated when omitted.

    Examples
    --------
    >>> MoleculeTrace(molecule="CCO", session_id="s", log_path=None).depth
    0
    """

    molecule: str
    session_id: str
    log_path: Path | None = None
    trace_id: str = field(default_factory=lambda: uuid.uuid4().hex)

    @property
    def depth(self) -> int:
        """Return the retrosynthesis depth of the node being expanded.

        Returns
        -------
        int
            Current node depth, ``0`` outside any :func:`node_depth` block.

        Examples
        --------
        >>> with molecule_trace("CCO", session_id="s") as trace:
        ...     with node_depth(3):
        ...         trace.depth
        3
        """
        return _depth_var.get()


def molecule_slug(smiles: str, *, limit: int = 60) -> str:
    """Build a filesystem- and session-id-safe name for a molecule.

    Mirrors :func:`deepretro.batch.slugify_molecule`, but lives here so that
    ``deepretro.batch`` (which imports this module) is not imported back.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.
    limit : int, optional
        Maximum length of the readable prefix.

    Returns
    -------
    str
        ASCII-safe prefix plus a short digest of the original SMILES.

    Examples
    --------
    >>> molecule_slug("CCO") == molecule_slug("CCO")
    True
    >>> "/" in molecule_slug("CC(=O)O")
    False
    """
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", smiles)[:limit]
    digest = hashlib.md5(smiles.encode("utf-8")).hexdigest()[:8]
    return f"{safe}_{digest}"


def elapsed_ms(started: float) -> float:
    """Return the milliseconds elapsed since a :func:`time.perf_counter` mark.

    Parameters
    ----------
    started : float
        Value previously returned by :func:`time.perf_counter`.

    Returns
    -------
    float
        Elapsed wall-clock time in milliseconds.

    Examples
    --------
    >>> import time
    >>> elapsed_ms(time.perf_counter()) >= 0.0
    True
    """
    return (time.perf_counter() - started) * 1000.0


def _default_session_id(molecule: str) -> str:
    """Build a readable, collision-resistant Langfuse session id.

    Examples
    --------
    >>> _default_session_id("CCO").startswith("autosolve_CCO_")
    True
    """
    stamp = time.strftime("%Y%m%d_%H%M%S")
    return f"{TRACE_NAME}_{molecule_slug(molecule, limit=40)}_{stamp}_{os.getpid()}"


@contextmanager
def molecule_trace(
    molecule: str,
    *,
    log_dir: str | Path | None = None,
    session_id: str | None = None,
) -> Iterator[MoleculeTrace]:
    """Activate the LLM trace context for one target molecule.

    **Nesting reuses the outer trace.** The batch runner opens a trace around
    each molecule so the log lands beside that molecule's pathway files, and
    :meth:`deepretro.algorithms.autosolve.AutoSolver.autosolve` opens one too so
    direct callers still get grouping. When a trace is already active this
    context manager yields it unchanged: ``log_dir`` and ``session_id`` are
    ignored, so the outer (batch) destination always wins and one run never
    splits into two Langfuse sessions.

    Parameters
    ----------
    molecule : str
        Target molecule SMILES.
    log_dir : str or Path or None, optional
        Directory that receives ``llm_calls.jsonl``. ``None`` (default)
        disables local logging and keeps only the Langfuse grouping.
    session_id : str or None, optional
        Explicit Langfuse session id. Generated as
        ``autosolve_<slug>_<YYYYmmdd_HHMMSS>_<pid>`` when omitted.

    Yields
    ------
    MoleculeTrace
        The active trace, which is the outer one when already nested.

    Examples
    --------
    >>> with molecule_trace("CCO", session_id="outer") as outer:
    ...     with molecule_trace("CCC", session_id="inner") as inner:
    ...         inner is outer
    True
    >>> current_trace() is None
    True
    """
    active = _trace_var.get()
    if active is not None:
        yield active
        return

    trace = MoleculeTrace(
        molecule=molecule,
        session_id=session_id or _default_session_id(molecule),
        log_path=None if log_dir is None else Path(log_dir) / LOG_FILENAME,
    )
    trace_token = _trace_var.set(trace)
    depth_token = _depth_var.set(0)
    node_token = _node_molecule_var.set(None)
    _open_langfuse_trace(trace)
    try:
        yield trace
    finally:
        _node_molecule_var.reset(node_token)
        _depth_var.reset(depth_token)
        _trace_var.reset(trace_token)
        _flush_langfuse()


def build_langfuse_client() -> Any | None:
    """Build (once) the Langfuse client used for tool events.

    Returns ``None`` when ``LANGFUSE_PUBLIC_KEY`` or ``LANGFUSE_SECRET_KEY`` is
    missing, when the ``langfuse`` package is unavailable, or when the client
    cannot be constructed. The result is cached for the process.

    Returns
    -------
    Any or None
        A ``langfuse.Langfuse`` instance, or ``None`` when disabled.

    Examples
    --------
    >>> import os
    >>> os.environ.pop("LANGFUSE_PUBLIC_KEY", None) and None
    >>> build_langfuse_client() is None  # doctest: +SKIP
    True
    """
    global _client_cache, _client_checked
    if _client_checked:
        return _client_cache
    _client_checked = True
    if not (os.getenv("LANGFUSE_PUBLIC_KEY") and os.getenv("LANGFUSE_SECRET_KEY")):
        return None
    try:
        from langfuse import Langfuse

        _client_cache = Langfuse()
    except Exception as exc:  # pragma: no cover - depends on the environment
        logger.warning("llm_trace.langfuse_client_unavailable", error=str(exc))
        _client_cache = None
    return _client_cache


def _langfuse_client() -> Any | None:
    """Indirection over :func:`build_langfuse_client` so tests can swap it.

    Examples
    --------
    >>> _langfuse_client() is build_langfuse_client()
    True
    """
    return build_langfuse_client()


def _open_langfuse_trace(trace: MoleculeTrace) -> None:
    """Create the Langfuse trace shell so tool events can attach to it.

    LiteLLM upserts the same trace id on every generation, so the shell only
    needs to exist; failures are logged once and ignored.

    Examples
    --------
    >>> _open_langfuse_trace(MoleculeTrace(molecule="CCO", session_id="s"))
    """
    client = _langfuse_client()
    if client is None:
        return
    try:
        client.trace(
            id=trace.trace_id,
            name=TRACE_NAME,
            session_id=trace.session_id,
            metadata={"molecule": trace.molecule},
            tags=["deepretro", TRACE_NAME],
        )
    except Exception as exc:
        _warn_once("llm_trace.langfuse_trace_failed", str(exc))


def _flush_langfuse() -> None:
    """Flush queued Langfuse events when a molecule trace closes.

    Examples
    --------
    >>> _flush_langfuse()
    """
    client = _langfuse_client()
    if client is None:
        return
    try:
        client.flush()
    except Exception as exc:
        _warn_once("llm_trace.langfuse_flush_failed", str(exc))


def current_trace() -> MoleculeTrace | None:
    """Return the active molecule trace, if any.

    Returns
    -------
    MoleculeTrace or None
        The active trace, or ``None`` outside :func:`molecule_trace`.

    Examples
    --------
    >>> current_trace() is None
    True
    >>> with molecule_trace("CCO", session_id="s"):
    ...     current_trace().molecule
    'CCO'
    """
    return _trace_var.get()


def current_depth() -> int:
    """Return the retrosynthesis depth of the node currently being expanded.

    Returns
    -------
    int
        Node depth, ``0`` outside any :func:`node_depth` block.

    Examples
    --------
    >>> current_depth()
    0
    """
    return _depth_var.get()


def current_node_molecule() -> str | None:
    """Return the molecule of the node currently being expanded.

    Returns
    -------
    str or None
        Node SMILES, or ``None`` when no :func:`node_depth` block set one.

    Examples
    --------
    >>> current_node_molecule() is None
    True
    """
    return _node_molecule_var.get()


@contextmanager
def node_depth(depth: int, *, molecule: str | None = None) -> Iterator[int]:
    """Mark the recursion depth (and molecule) of the node being expanded.

    Wrapping a recursive LLM call in this block makes every record written by
    :func:`record_llm_call` carry where in the retrosynthesis tree it happened.

    Parameters
    ----------
    depth : int
        Recursion depth of the node.
    molecule : str or None, optional
        SMILES of the node being expanded, which differs from the trace's
        target molecule below the root.

    Yields
    ------
    int
        The depth that was set.

    Examples
    --------
    >>> with node_depth(2, molecule="CCC"):
    ...     (current_depth(), current_node_molecule())
    (2, 'CCC')
    >>> current_depth()
    0
    """
    depth_token = _depth_var.set(depth)
    node_token = _node_molecule_var.set(molecule)
    try:
        yield depth
    finally:
        _node_molecule_var.reset(node_token)
        _depth_var.reset(depth_token)


def langfuse_metadata(
    base: dict[str, Any] | None,
    *,
    stage: str,
) -> dict[str, Any]:
    """Merge caller metadata with the active trace's Langfuse grouping keys.

    Adds the keys LiteLLM's Langfuse integration reads (``trace_id``,
    ``session_id``, ``trace_name``, ``generation_name``, ``trace_metadata``,
    ``tags``) plus ``depth`` and ``node_molecule``, which LiteLLM forwards as
    per-generation metadata so each generation shows where in the tree it
    happened.

    Parameters
    ----------
    base : dict or None
        Caller metadata, for example ``{"task": "retrosynthesis"}``. Never
        modified in place.
    stage : str
        Which call site this is (``retrosynthesis``, ``retrosynthesis_agent``,
        ``metadata``). Used as the Langfuse generation name.

    Returns
    -------
    dict
        A new mapping. Without an active trace it is ``base`` plus ``stage``.

    Examples
    --------
    >>> langfuse_metadata({"task": "metadata"}, stage="metadata")
    {'task': 'metadata', 'stage': 'metadata'}
    >>> with molecule_trace("CCO", session_id="s"):
    ...     langfuse_metadata(None, stage="metadata")["trace_name"]
    'autosolve'
    """
    metadata: dict[str, Any] = dict(base or {})
    metadata["stage"] = stage

    trace = _trace_var.get()
    if trace is None:
        return metadata

    node_molecule = _node_molecule_var.get()
    depth = _depth_var.get()
    metadata.update(
        {
            "trace_id": trace.trace_id,
            "session_id": trace.session_id,
            "trace_name": TRACE_NAME,
            "generation_name": stage,
            "trace_metadata": {
                "molecule": trace.molecule,
                "node_molecule": node_molecule,
                "depth": depth,
            },
            "tags": ["deepretro", TRACE_NAME, stage],
            "depth": depth,
            "node_molecule": node_molecule,
        }
    )
    return metadata


def serialize_tool_calls(tool_calls: Any) -> Any:
    """Turn provider tool-call objects into plain JSON-friendly dicts.

    Parameters
    ----------
    tool_calls : Any
        ``None``, or a sequence of dicts / Pydantic-like objects exposing
        ``model_dump()``.

    Returns
    -------
    Any
        ``None``, or a list of dicts (``str`` for anything unrecognised).

    Examples
    --------
    >>> serialize_tool_calls(None) is None
    True
    >>> serialize_tool_calls([{"id": "c1"}])
    [{'id': 'c1'}]
    """
    if tool_calls is None:
        return None
    serialized: list[Any] = []
    for call in tool_calls:
        if isinstance(call, dict):
            serialized.append(call)
            continue
        dump = getattr(call, "model_dump", None)
        if callable(dump):
            try:
                serialized.append(dump())
                continue
            except Exception:  # pragma: no cover - provider serializers may fail
                pass
        serialized.append(str(call))
    return serialized


def sanitize_messages(messages: Sequence[Any] | None) -> list[Any] | None:
    """Copy a conversation for logging, without opaque provider blobs.

    Anthropic returns signed ``thinking_blocks`` that must be replayed on the
    wire but are several kilobytes of base64 per turn with no readable
    content. They and ``provider_specific_fields`` are dropped, an empty
    ``reasoning_content`` is dropped, and tool calls are converted to plain
    dicts. The caller's messages are never modified.

    Parameters
    ----------
    messages : sequence or None
        Conversation as sent to the provider.

    Returns
    -------
    list or None
        Cleaned copies, or ``None`` when ``messages`` is ``None``.

    Examples
    --------
    >>> sanitize_messages([{"role": "assistant", "content": "x",
    ...                     "thinking_blocks": [{"signature": "..."}]}])
    [{'role': 'assistant', 'content': 'x'}]
    """
    if messages is None:
        return None
    cleaned: list[Any] = []
    for message in messages:
        dump = getattr(message, "model_dump", None)
        if not isinstance(message, dict) and callable(dump):
            try:
                message = dump()
            except Exception:  # pragma: no cover - provider serializers may fail
                cleaned.append(str(message))
                continue
        if not isinstance(message, dict):
            cleaned.append(str(message))
            continue
        copy = {k: v for k, v in message.items() if k not in _OPAQUE_MESSAGE_KEYS}
        if not copy.get("reasoning_content"):
            copy.pop("reasoning_content", None)
        if copy.get("tool_calls") is not None:
            copy["tool_calls"] = serialize_tool_calls(copy["tool_calls"])
        cleaned.append(copy)
    return cleaned


def _response_message(response: Any) -> Any:
    """Pull the assistant message out of a LiteLLM response, if present.

    Examples
    --------
    >>> _response_message({"role": "assistant"})
    {'role': 'assistant'}
    """
    choices = getattr(response, "choices", None)
    if choices:
        return getattr(choices[0], "message", choices[0])
    return response


def _serialize_response(response: Any) -> Any:
    """Reduce a LiteLLM response to its content string or serialized message.

    Examples
    --------
    >>> _serialize_response(None) is None
    True
    >>> _serialize_response("done")
    'done'
    """
    if response is None or isinstance(response, str):
        return response

    message = _response_message(response)
    if isinstance(message, dict):
        content = message.get("content")
        return content if isinstance(content, str) and content else message

    content = getattr(message, "content", None)
    if isinstance(content, str) and content:
        return content

    dump = getattr(message, "model_dump", None)
    if callable(dump):
        try:
            return dump()
        except Exception:  # pragma: no cover - provider serializers may fail
            return str(message)
    return str(message)


def _extract_usage(response: Any) -> dict[str, int] | None:
    """Read prompt/completion/total token counts off a LiteLLM response.

    Examples
    --------
    >>> _extract_usage(None) is None
    True
    """
    usage = getattr(response, "usage", None)
    if usage is None:
        return None
    counts: dict[str, int] = {}
    for key in ("prompt_tokens", "completion_tokens", "total_tokens"):
        value = usage.get(key) if isinstance(usage, dict) else getattr(usage, key, None)
        if isinstance(value, int):
            counts[key] = value
    return counts or None


def record_llm_call(
    *,
    stage: str,
    model: str,
    messages: Sequence[Any] | None,
    response: Any,
    error: str | None = None,
    latency_ms: float,
    usage: dict[str, Any] | None = None,
    tool_calls: Any = None,
    iteration: int | None = None,
    node_molecule: str | None = None,
) -> None:
    """Append one LLM call to the active trace's local JSONL log.

    A no-op when no trace is active or the trace has no log path. Writing is
    best-effort: serialization falls back to ``str`` and filesystem errors are
    logged once and swallowed, so observability can never break a run.

    Parameters
    ----------
    stage : str
        Call site (``retrosynthesis``, ``retrosynthesis_agent``, ``metadata``).
    model : str
        LiteLLM model identifier used for the call.
    messages : sequence or None
        Conversation sent to the provider.
    response : Any
        LiteLLM response, assistant message, or content string. ``None`` on
        failure.
    error : str or None, optional
        Error text when the call failed.
    latency_ms : float
        Wall-clock duration of the provider call.
    usage : dict or None, optional
        Token usage override; read off ``response`` when omitted.
    tool_calls : Any, optional
        Tool calls requested by the assistant turn.
    iteration : int or None, optional
        Agent loop iteration or provider retry attempt, 1-based.
    node_molecule : str or None, optional
        Molecule of the node being expanded. Defaults to the value set by
        :func:`node_depth`.

    Examples
    --------
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     with molecule_trace("CCO", log_dir=tmp, session_id="s") as trace:
    ...         record_llm_call(
    ...             stage="retrosynthesis",
    ...             model="openai/gpt-4o-mini",
    ...             messages=[{"role": "user", "content": "hi"}],
    ...             response="done",
    ...             latency_ms=1.0,
    ...         )
    ...     len(trace.log_path.read_text().splitlines())
    1
    """
    trace = _trace_var.get()
    if trace is None or trace.log_path is None:
        return

    # Serialization touches provider objects, so keep it inside the guard: an
    # exotic response must never abort the LLM call it is describing.
    try:
        record = _envelope(trace, "llm_call", stage, iteration, node_molecule)
        record.update(
            {
                "model": model,
                "messages": sanitize_messages(messages),
                "response": _serialize_response(response),
                "tool_calls": serialize_tool_calls(tool_calls),
                "usage": usage if usage is not None else _extract_usage(response),
                "latency_ms": latency_ms,
                "error": error,
            }
        )
    except Exception as exc:  # pragma: no cover - default=str makes this rare
        _warn_once("llm_trace.serialize_failed", str(exc))
        return
    _write_record(trace, record)


def record_tool_results(
    *,
    stage: str,
    iteration: int | None,
    results: Sequence[dict[str, Any]],
) -> None:
    """Log the tools the agent ran after one model turn.

    Writes one ``kind="tool_results"`` line to the local log and, when a
    Langfuse client is configured, one Langfuse *event* per tool inside the
    molecule's trace (``name="tool:<tool name>"``, ``input`` = arguments,
    ``output`` = tool result) so the outputs render online beside the model
    turns. A no-op without an active trace. Never raises.

    Parameters
    ----------
    stage : str
        Call site, normally ``retrosynthesis_agent``.
    iteration : int or None
        1-based agent turn whose tool calls these results answer.
    results : sequence of dict
        One mapping per executed tool with keys ``tool_call_id``, ``name``,
        ``arguments`` and ``output``.

    Examples
    --------
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     with molecule_trace("CCO", log_dir=tmp, session_id="s") as trace:
    ...         record_tool_results(
    ...             stage="retrosynthesis_agent",
    ...             iteration=1,
    ...             results=[{"tool_call_id": "c", "name": "validate_smiles",
    ...                       "arguments": {"smiles": "CCO"},
    ...                       "output": {"valid": True}}],
    ...         )
    ...     json.loads(trace.log_path.read_text())["kind"]
    'tool_results'
    """
    trace = _trace_var.get()
    if trace is None:
        return
    results = list(results)

    if trace.log_path is not None:
        try:
            record = _envelope(trace, "tool_results", stage, iteration, None)
            record["tool_results"] = results
        except Exception as exc:  # pragma: no cover - plain dicts; defensive
            _warn_once("llm_trace.serialize_failed", str(exc))
        else:
            _write_record(trace, record)

    client = _langfuse_client()
    if client is None:
        return
    depth = _depth_var.get()
    node_molecule = _node_molecule_var.get()
    for result in results:
        try:
            client.event(
                trace_id=trace.trace_id,
                name=f"tool:{result.get('name', 'unknown')}",
                input=result.get("arguments"),
                output=result.get("output"),
                metadata={
                    "stage": stage,
                    "iteration": iteration,
                    "depth": depth,
                    "node_molecule": node_molecule,
                    "tool_call_id": result.get("tool_call_id"),
                },
            )
        except Exception as exc:
            _warn_once("llm_trace.langfuse_event_failed", str(exc))


def _envelope(
    trace: MoleculeTrace,
    kind: str,
    stage: str,
    iteration: int | None,
    node_molecule: str | None,
) -> dict[str, Any]:
    """Build the fields every log line shares.

    Examples
    --------
    >>> t = MoleculeTrace(molecule="CCO", session_id="s")
    >>> _envelope(t, "llm_call", "retrosynthesis", 1, None)["kind"]
    'llm_call'
    """
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "kind": kind,
        "session_id": trace.session_id,
        "trace_id": trace.trace_id,
        "target": trace.molecule,
        "node_molecule": node_molecule or _node_molecule_var.get(),
        "depth": _depth_var.get(),
        "stage": stage,
        "iteration": iteration,
    }


def _write_record(trace: MoleculeTrace, record: dict[str, Any]) -> None:
    """Append one JSON line to the trace's log; filesystem errors warn once.

    Examples
    --------
    >>> _write_record(MoleculeTrace(molecule="CCO", session_id="s"), {})
    """
    if trace.log_path is None:
        return
    try:
        line = json.dumps(record, default=str)
        trace.log_path.parent.mkdir(parents=True, exist_ok=True)
        with trace.log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
    except (OSError, TypeError, ValueError) as exc:
        _warn_once("llm_trace.write_failed", str(exc))


def _warn_once(event: str, error: str) -> None:
    """Log a trace-logging failure once per process.

    Examples
    --------
    >>> _warn_once("llm_trace.write_failed", "disk full")  # doctest: +SKIP
    """
    global _write_failure_logged
    if _write_failure_logged:
        return
    _write_failure_logged = True
    logger.warning(event, error=error)
