"""Tests for the per-molecule LLM trace context and local JSONL call log."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from deepretro.utils import llm_trace
from deepretro.utils.llm_trace import (
    LOG_FILENAME,
    MoleculeTrace,
    current_trace,
    elapsed_ms,
    langfuse_metadata,
    molecule_trace,
    node_depth,
    record_llm_call,
)


class FakeMessage:
    """Assistant message stand-in exposing content and a Pydantic serializer."""

    def __init__(self, content: str | None = "ok") -> None:
        self.content = content
        self.tool_calls: list[dict[str, Any]] | None = None

    def model_dump(self) -> dict[str, Any]:
        return {"role": "assistant", "content": self.content}


class FakeChoice:
    """LiteLLM choice stand-in."""

    def __init__(self, message: FakeMessage) -> None:
        self.message = message


class FakeUsage:
    """LiteLLM usage stand-in."""

    prompt_tokens = 11
    completion_tokens = 7
    total_tokens = 18


class FakeResponse:
    """Minimal LiteLLM ``ModelResponse`` stand-in."""

    def __init__(self, content: str | None = "ok", usage: Any = None) -> None:
        self.choices = [FakeChoice(FakeMessage(content))]
        self.usage = usage if usage is not None else FakeUsage()


def read_records(log_dir: Path) -> list[dict[str, Any]]:
    """Read every JSON line written to ``log_dir/llm_calls.jsonl``."""
    path = log_dir / LOG_FILENAME
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


# ----------------------------------------------------------------------
# Context manager
# ----------------------------------------------------------------------


def test_no_trace_is_active_by_default() -> None:
    """Outside any context manager there is no active trace."""
    assert current_trace() is None


def test_molecule_trace_sets_and_clears(tmp_path: Path) -> None:
    """The context manager activates a trace and restores the previous state."""
    with molecule_trace("CCO", log_dir=tmp_path) as trace:
        assert current_trace() is trace
        assert trace.molecule == "CCO"
        assert trace.log_path == tmp_path / LOG_FILENAME
    assert current_trace() is None


def test_molecule_trace_clears_on_exception(tmp_path: Path) -> None:
    """A raising body still resets the trace context."""
    with pytest.raises(RuntimeError), molecule_trace("CCO", log_dir=tmp_path):
        raise RuntimeError("boom")
    assert current_trace() is None


def test_session_id_is_generated_and_slugged() -> None:
    """A generated session id carries the molecule slug and the process id."""
    import os

    with molecule_trace("CC(=O)O") as trace:
        assert trace.session_id.startswith("autosolve_")
        assert "/" not in trace.session_id
        assert trace.session_id.endswith(f"_{os.getpid()}")


def test_explicit_session_id_is_used() -> None:
    """An explicit session id overrides generation."""
    with molecule_trace("CCO", session_id="fixed") as trace:
        assert trace.session_id == "fixed"


def test_log_path_is_none_without_log_dir() -> None:
    """Without a log directory the trace only groups Langfuse generations."""
    with molecule_trace("CCO") as trace:
        assert trace.log_path is None


def test_nesting_reuses_the_outer_trace(tmp_path: Path) -> None:
    """An inner context manager must not replace an active outer trace."""
    inner_dir = tmp_path / "inner"
    with molecule_trace("CCO", log_dir=tmp_path, session_id="outer") as outer:
        with molecule_trace("CCC", log_dir=inner_dir, session_id="inner") as inner:
            assert inner is outer
            assert inner.session_id == "outer"
            assert inner.molecule == "CCO"
            assert inner.log_path == tmp_path / LOG_FILENAME
        assert current_trace() is outer
    assert current_trace() is None


# ----------------------------------------------------------------------
# Depth context
# ----------------------------------------------------------------------


def test_depth_defaults_to_zero_and_nests() -> None:
    """``node_depth`` sets the current depth and restores it on exit."""
    with molecule_trace("CCO") as trace:
        assert trace.depth == 0
        with node_depth(2):
            assert trace.depth == 2
            with node_depth(3):
                assert trace.depth == 3
            assert trace.depth == 2
        assert trace.depth == 0


def test_node_depth_carries_the_node_molecule() -> None:
    """``node_depth`` may also record which molecule the node is expanding."""
    with molecule_trace("CCO"):
        with node_depth(1, molecule="CCC"):
            assert llm_trace.current_node_molecule() == "CCC"
        assert llm_trace.current_node_molecule() is None


# ----------------------------------------------------------------------
# langfuse_metadata
# ----------------------------------------------------------------------


def test_langfuse_metadata_is_a_no_op_without_a_trace() -> None:
    """Without an active trace only the stage is added."""
    metadata = langfuse_metadata({"task": "retrosynthesis"}, stage="retrosynthesis")
    assert metadata == {"task": "retrosynthesis", "stage": "retrosynthesis"}


def test_langfuse_metadata_accepts_none_base() -> None:
    """A ``None`` base is treated as an empty mapping."""
    assert langfuse_metadata(None, stage="metadata") == {"stage": "metadata"}


def test_langfuse_metadata_merges_trace_fields() -> None:
    """With an active trace the LiteLLM/Langfuse grouping keys are added."""
    with molecule_trace("CCO", session_id="sess-1"):
        with node_depth(2, molecule="CCC"):
            metadata = langfuse_metadata({"task": "x"}, stage="retrosynthesis_agent")

    assert metadata["task"] == "x"
    assert metadata["session_id"] == "sess-1"
    assert metadata["trace_name"] == "autosolve"
    assert metadata["generation_name"] == "retrosynthesis_agent"
    assert metadata["trace_metadata"] == {
        "molecule": "CCO",
        "node_molecule": "CCC",
        "depth": 2,
    }
    assert "retrosynthesis_agent" in metadata["tags"]


def test_langfuse_metadata_does_not_mutate_the_base() -> None:
    """The caller's mapping is never modified in place."""
    base = {"task": "x"}
    with molecule_trace("CCO", session_id="sess-1"):
        langfuse_metadata(base, stage="metadata")
    assert base == {"task": "x"}


# ----------------------------------------------------------------------
# record_llm_call
# ----------------------------------------------------------------------


def test_record_llm_call_is_a_no_op_without_a_trace(tmp_path: Path) -> None:
    """No trace means no file is created."""
    record_llm_call(
        stage="retrosynthesis",
        model="openai/gpt-4o-mini",
        messages=[{"role": "user", "content": "hi"}],
        response=FakeResponse(),
        latency_ms=1.0,
    )
    assert not (tmp_path / LOG_FILENAME).exists()


def test_record_llm_call_is_a_no_op_without_a_log_path() -> None:
    """A trace without a log directory records nothing locally."""
    with molecule_trace("CCO") as trace:
        assert trace.log_path is None
        record_llm_call(
            stage="retrosynthesis",
            model="openai/gpt-4o-mini",
            messages=[],
            response=FakeResponse(),
            latency_ms=1.0,
        )


def test_record_llm_call_writes_expected_fields(tmp_path: Path) -> None:
    """One call appends one JSON line carrying every documented field."""
    with molecule_trace("CCO", log_dir=tmp_path, session_id="sess-1"):
        with node_depth(1, molecule="CCC"):
            record_llm_call(
                stage="retrosynthesis",
                model="openai/gpt-4o-mini",
                messages=[{"role": "user", "content": "hi"}],
                response=FakeResponse("answer"),
                latency_ms=12.5,
                iteration=2,
            )

    records = read_records(tmp_path)
    assert len(records) == 1
    record = records[0]
    assert record["session_id"] == "sess-1"
    assert record["target"] == "CCO"
    assert record["node_molecule"] == "CCC"
    assert record["depth"] == 1
    assert record["stage"] == "retrosynthesis"
    assert record["iteration"] == 2
    assert record["model"] == "openai/gpt-4o-mini"
    assert record["messages"] == [{"role": "user", "content": "hi"}]
    assert record["response"] == "answer"
    assert record["latency_ms"] == 12.5
    assert record["error"] is None
    assert record["usage"] == {
        "prompt_tokens": 11,
        "completion_tokens": 7,
        "total_tokens": 18,
    }
    assert record["timestamp"].endswith("+00:00")


def test_record_llm_call_appends_one_line_per_call(tmp_path: Path) -> None:
    """Records accumulate in the same file across calls."""
    with molecule_trace("CCO", log_dir=tmp_path):
        for index in range(3):
            record_llm_call(
                stage="retrosynthesis_agent",
                model="m",
                messages=[],
                response=FakeResponse(f"r{index}"),
                latency_ms=1.0,
                iteration=index,
            )
    assert [record["response"] for record in read_records(tmp_path)] == [
        "r0",
        "r1",
        "r2",
    ]


def test_record_llm_call_serializes_a_message_without_content(tmp_path: Path) -> None:
    """A tool-call-only message is stored as its serialized mapping."""
    with molecule_trace("CCO", log_dir=tmp_path):
        record_llm_call(
            stage="retrosynthesis_agent",
            model="m",
            messages=[],
            response=FakeResponse(None),
            latency_ms=1.0,
            tool_calls=[{"id": "call_1", "function": {"name": "validate_smiles"}}],
        )
    record = read_records(tmp_path)[0]
    assert record["response"] == {"role": "assistant", "content": None}
    assert record["tool_calls"][0]["id"] == "call_1"


def test_record_llm_call_records_errors(tmp_path: Path) -> None:
    """A failed call is recorded with its error and a null response."""
    with molecule_trace("CCO", log_dir=tmp_path):
        record_llm_call(
            stage="metadata",
            model="m",
            messages=[],
            response=None,
            error="rate limited",
            latency_ms=3.0,
        )
    record = read_records(tmp_path)[0]
    assert record["error"] == "rate limited"
    assert record["response"] is None
    assert record["usage"] is None


def test_record_llm_call_serializes_unusual_objects(tmp_path: Path) -> None:
    """Non-JSON-serializable payloads are stringified instead of raising."""

    class Odd:
        def __repr__(self) -> str:
            return "<odd>"

    with molecule_trace("CCO", log_dir=tmp_path):
        record_llm_call(
            stage="retrosynthesis",
            model="m",
            messages=[{"role": "user", "content": Odd()}],
            response="text",
            latency_ms=1.0,
        )
    record = read_records(tmp_path)[0]
    assert record["messages"][0]["content"] == "<odd>"


def test_record_llm_call_swallows_os_errors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A filesystem failure must never propagate into the solver."""

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise OSError("disk full")

    with molecule_trace("CCO", log_dir=tmp_path):
        monkeypatch.setattr(Path, "open", boom)
        record_llm_call(
            stage="retrosynthesis",
            model="m",
            messages=[],
            response="text",
            latency_ms=1.0,
        )


def test_record_llm_call_prefers_explicit_usage(tmp_path: Path) -> None:
    """An explicitly supplied usage mapping wins over the response usage."""
    with molecule_trace("CCO", log_dir=tmp_path):
        record_llm_call(
            stage="retrosynthesis",
            model="m",
            messages=[],
            response=FakeResponse(),
            latency_ms=1.0,
            usage={"total_tokens": 1},
        )
    assert read_records(tmp_path)[0]["usage"] == {"total_tokens": 1}


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def test_elapsed_ms_is_non_negative() -> None:
    """``elapsed_ms`` converts a perf counter start into milliseconds."""
    import time

    assert elapsed_ms(time.perf_counter()) >= 0.0


def test_molecule_trace_dataclass_is_constructible() -> None:
    """``MoleculeTrace`` is a plain value object."""
    trace = MoleculeTrace(molecule="CCO", session_id="s", log_path=None)
    assert trace.molecule == "CCO"
    assert trace.depth == 0


# ---------------------------------------------------------------------------
# Record slimming, tool results, Langfuse events
# ---------------------------------------------------------------------------


class FakeToolCall:
    """Pydantic-like tool call object as LiteLLM returns it."""

    def __init__(self, call_id: str, name: str, arguments: str) -> None:
        self.id = call_id
        self.function = {"name": name, "arguments": arguments}

    def model_dump(self) -> dict[str, Any]:
        return {"id": self.id, "type": "function", "function": self.function}


class FakeLangfuse:
    """Recorder standing in for the Langfuse v2 client."""

    def __init__(self) -> None:
        self.traces: list[dict[str, Any]] = []
        self.events: list[dict[str, Any]] = []
        self.flushed = 0

    def trace(self, **kwargs: Any) -> None:
        self.traces.append(kwargs)

    def event(self, **kwargs: Any) -> None:
        self.events.append(kwargs)

    def flush(self) -> None:
        self.flushed += 1


def _read(tmp_path: Path) -> list[dict[str, Any]]:
    path = tmp_path / LOG_FILENAME
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def test_llm_call_records_carry_kind_and_drop_thinking_signatures(
    tmp_path: Path,
) -> None:
    messages = [
        {"role": "user", "content": "hi"},
        {
            "role": "assistant",
            "content": None,
            "tool_calls": [FakeToolCall("c1", "validate_smiles", '{"smiles": "CCO"}')],
            "reasoning_content": "",
            "thinking_blocks": [{"type": "thinking", "signature": "x" * 5000}],
            "provider_specific_fields": {"thinking_blocks": [{"signature": "y"}]},
        },
    ]
    with molecule_trace("CCO", log_dir=tmp_path, session_id="s"):
        record_llm_call(
            stage="retrosynthesis_agent",
            model="m",
            messages=messages,
            response="done",
            latency_ms=1.0,
            tool_calls=[FakeToolCall("c2", "check_hallucination", "{}")],
        )
    (record,) = _read(tmp_path)
    assert record["kind"] == "llm_call"
    assistant = record["messages"][1]
    assert "thinking_blocks" not in assistant
    assert "provider_specific_fields" not in assistant
    assert "reasoning_content" not in assistant  # empty string dropped
    assert assistant["tool_calls"][0] == {
        "id": "c1",
        "type": "function",
        "function": {"name": "validate_smiles", "arguments": '{"smiles": "CCO"}'},
    }
    assert record["tool_calls"][0]["function"]["name"] == "check_hallucination"
    assert len(json.dumps(record)) < 1500
    # The caller's message list is not mutated.
    assert "thinking_blocks" in messages[1]


def test_non_empty_reasoning_content_is_kept() -> None:
    cleaned = llm_trace.sanitize_messages(
        [{"role": "assistant", "content": "x", "reasoning_content": "because"}]
    )
    assert cleaned[0]["reasoning_content"] == "because"


def test_record_tool_results_writes_a_tool_results_line(tmp_path: Path) -> None:
    with molecule_trace("CCO", log_dir=tmp_path, session_id="s"):
        with node_depth(1, molecule="CC=O"):
            llm_trace.record_tool_results(
                stage="retrosynthesis_agent",
                iteration=2,
                results=[
                    {
                        "tool_call_id": "c1",
                        "name": "validate_smiles",
                        "arguments": {"smiles": "CCO"},
                        "output": {"valid": True},
                    }
                ],
            )
    (record,) = _read(tmp_path)
    assert record["kind"] == "tool_results"
    assert record["stage"] == "retrosynthesis_agent"
    assert record["iteration"] == 2
    assert record["depth"] == 1
    assert record["node_molecule"] == "CC=O"
    assert record["tool_results"][0]["name"] == "validate_smiles"
    assert record["tool_results"][0]["output"] == {"valid": True}


def test_record_tool_results_is_a_noop_without_a_trace(tmp_path: Path) -> None:
    llm_trace.record_tool_results(
        stage="retrosynthesis_agent", iteration=1, results=[{"name": "x"}]
    )
    assert not (tmp_path / LOG_FILENAME).exists()


def test_trace_has_a_stable_langfuse_trace_id_used_by_metadata() -> None:
    with molecule_trace("CCO", session_id="s") as trace:
        assert trace.trace_id
        with node_depth(2, molecule="CC=O"):
            meta = langfuse_metadata({"task": "t"}, stage="retrosynthesis_agent")
    assert meta["trace_id"] == trace.trace_id
    # Per-generation keys so Langfuse shows them on each generation online.
    assert meta["depth"] == 2
    assert meta["node_molecule"] == "CC=O"
    assert meta["trace_metadata"] == {
        "molecule": "CCO",
        "node_molecule": "CC=O",
        "depth": 2,
    }


def test_tool_results_become_langfuse_events_in_the_same_trace(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = FakeLangfuse()
    monkeypatch.setattr(llm_trace, "_langfuse_client", lambda: fake)

    with molecule_trace("CCO", log_dir=tmp_path, session_id="s") as trace:
        with node_depth(1, molecule="CC=O"):
            llm_trace.record_tool_results(
                stage="retrosynthesis_agent",
                iteration=3,
                results=[
                    {
                        "tool_call_id": "c1",
                        "name": "validate_smiles",
                        "arguments": {"smiles": "CCO"},
                        "output": {"valid": True},
                    },
                    {
                        "tool_call_id": "c2",
                        "name": "check_hallucination",
                        "arguments": {"product": "CC=O"},
                        "output": {"score": 90},
                    },
                ],
            )

    # The trace shell is created up front so events can attach to it.
    assert fake.traces and fake.traces[0]["id"] == trace.trace_id
    assert fake.traces[0]["session_id"] == "s"
    assert fake.traces[0]["name"] == llm_trace.TRACE_NAME
    assert [e["name"] for e in fake.events] == [
        "tool:validate_smiles",
        "tool:check_hallucination",
    ]
    event = fake.events[0]
    assert event["trace_id"] == trace.trace_id
    assert event["input"] == {"smiles": "CCO"}
    assert event["output"] == {"valid": True}
    assert event["metadata"]["depth"] == 1
    assert event["metadata"]["iteration"] == 3
    assert event["metadata"]["tool_call_id"] == "c1"
    assert fake.flushed == 1  # flushed once when the trace closes


def test_langfuse_client_failures_never_propagate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class Broken:
        def trace(self, **kwargs: Any) -> None:
            raise RuntimeError("boom")

        def event(self, **kwargs: Any) -> None:
            raise RuntimeError("boom")

        def flush(self) -> None:
            raise RuntimeError("boom")

    monkeypatch.setattr(llm_trace, "_langfuse_client", lambda: Broken())
    with molecule_trace("CCO", log_dir=tmp_path, session_id="s"):
        llm_trace.record_tool_results(
            stage="retrosynthesis_agent",
            iteration=1,
            results=[{"tool_call_id": "c", "name": "n", "arguments": {}, "output": 1}],
        )
    assert len(_read(tmp_path)) == 1


def test_no_langfuse_client_without_keys(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("LANGFUSE_PUBLIC_KEY", raising=False)
    monkeypatch.delenv("LANGFUSE_SECRET_KEY", raising=False)
    monkeypatch.setattr(llm_trace, "_client_cache", None)
    assert llm_trace.build_langfuse_client() is None
