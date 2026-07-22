"""Tests for the single-step tool-calling agent loop."""

from __future__ import annotations

import json
from typing import Any

import pytest

import deepretro.agents.loop as agent_loop
from deepretro.agents.loop import agentic_orchestrator, agentic_single_step

MODEL = "openai/gpt-4o-mini"
FINAL_ANSWER = (
    '<json>{"data": [["CCO"]], "explanation": ["reduce"], '
    '"confidence_scores": [0.8]}</json>'
)


def final_turn(content: str = FINAL_ANSWER) -> dict[str, Any]:
    """An assistant turn with no tool calls (a final answer)."""
    return {"role": "assistant", "content": content}


def tool_turn(
    name: str = "validate_smiles",
    arguments: dict[str, Any] | None = None,
    call_id: str = "call_1",
) -> dict[str, Any]:
    """An assistant turn requesting one tool call."""
    return {
        "role": "assistant",
        "content": None,
        "tool_calls": [
            {
                "id": call_id,
                "type": "function",
                "function": {
                    "name": name,
                    "arguments": json.dumps(arguments or {"smiles": "CCO"}),
                },
            }
        ],
    }


class ScriptedModel:
    """Injectable model that replays scripted turns and records inputs."""

    def __init__(self, turns: list[dict[str, Any]]) -> None:
        self.turns = list(turns)
        self.seen: list[list[dict[str, Any]]] = []

    def __call__(self, messages: list[dict[str, Any]]) -> dict[str, Any]:
        self.seen.append([dict(message) for message in messages])
        return self.turns[len(self.seen) - 1]


class ProviderMessage:
    """Provider-message stand-in exposing a Pydantic-style serializer."""

    def __init__(self, payload: dict[str, Any]) -> None:
        self.payload = payload

    def model_dump(self) -> dict[str, Any]:
        """Return the provider payload."""
        return self.payload


def test_returns_parsed_pathways_from_final_answer() -> None:
    """A final answer is parsed into the pipeline's (pathways, expl, conf) shape."""
    model = ScriptedModel([final_turn()])
    pathways, explanations, confidence = agentic_single_step(
        "CC=O", MODEL, llm_runner=model
    )
    assert pathways == [["CCO"]]
    assert explanations == ["reduce"]
    assert confidence == [0.8]


def test_executes_tool_then_returns_final_answer() -> None:
    """A tool call is executed and its result fed back before the final answer."""
    model = ScriptedModel([tool_turn(), final_turn()])
    pathways, _explanations, _confidence = agentic_single_step(
        "CC=O", MODEL, llm_runner=model
    )
    assert pathways == [["CCO"]]
    second_call_messages = model.seen[1]
    assert any(message.get("role") == "tool" for message in second_call_messages)


def test_replays_provider_message_with_reasoning_fields() -> None:
    """A serialized provider turn is replayed with opaque reasoning intact."""
    first_payload = tool_turn()
    first_payload["reasoning_content"] = {"signature": "signed-token"}
    seen_second_turn: list[dict[str, Any]] = []
    calls = 0

    def provider_model(messages: list[dict[str, Any]]) -> Any:
        nonlocal calls
        calls += 1
        if calls == 1:
            return ProviderMessage(first_payload)
        first_payload["reasoning_content"]["signature"] = "changed"
        seen_second_turn.extend(messages)
        return ProviderMessage(final_turn())

    result = agentic_single_step("CC=O", MODEL, llm_runner=provider_model)

    assistant_turn = next(
        message for message in seen_second_turn if message["role"] == "assistant"
    )
    assert assistant_turn["reasoning_content"] == {"signature": "signed-token"}
    assert result == ([["CCO"]], ["reduce"], [0.8])


def test_thinking_is_enabled_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    """The default model-call factory receives thinking enabled."""
    captured: dict[str, Any] = {}

    def fake_model_call_factory(
        model: str,
        tools: list[dict[str, Any]],
        enable_thinking: bool,
        max_output_tokens: int | None,
    ) -> Any:
        captured["enable_thinking"] = enable_thinking
        return lambda messages: final_turn()

    monkeypatch.setattr(agent_loop, "_make_default_model_call", fake_model_call_factory)

    agentic_single_step("CC=O", MODEL)

    assert captured["enable_thinking"] is True


def test_tool_result_references_call_id() -> None:
    """The tool result message references the originating tool_call id."""
    model = ScriptedModel([tool_turn(call_id="abc123"), final_turn()])
    agentic_single_step("CC=O", MODEL, llm_runner=model)
    tool_messages = [m for m in model.seen[1] if m.get("role") == "tool"]
    assert tool_messages[0]["tool_call_id"] == "abc123"


def test_hits_max_iterations_returns_empty() -> None:
    """Looping without a final answer returns empty results after the cap."""

    def always_tool(messages: list[dict[str, Any]]) -> dict[str, Any]:
        return tool_turn()

    result = agentic_single_step(
        "CC=O", MODEL, llm_runner=always_tool, max_iterations=3
    )
    assert result == ([], [], [])


def test_unparseable_final_answer_returns_empty() -> None:
    """A final answer without a valid JSON payload yields empty results."""
    model = ScriptedModel([final_turn(content="I could not find a route.")])
    assert agentic_single_step("CC=O", MODEL, llm_runner=model) == ([], [], [])


def test_event_sink_records_refusal() -> None:
    """A refusal final message is recorded as a 'refusal' event."""
    model = ScriptedModel([final_turn(content="I can't help with that request.")])
    sink: list[dict[str, Any]] = []
    result = agentic_single_step("CC=O", MODEL, llm_runner=model, event_sink=sink)
    assert result == ([], [], [])
    assert len(sink) == 1
    assert sink[0]["kind"] == "refusal"
    assert sink[0]["molecule"] == "CC=O"


def test_event_sink_records_no_parseable_answer() -> None:
    """A non-refusal, unparseable final message is a 'no_parseable_answer' event."""
    model = ScriptedModel([final_turn(content="Here is a route but no JSON payload.")])
    sink: list[dict[str, Any]] = []
    agentic_single_step("CC=O", MODEL, llm_runner=model, event_sink=sink)
    assert len(sink) == 1
    assert sink[0]["kind"] == "no_parseable_answer"


def test_event_sink_records_no_final_answer_on_max_iterations() -> None:
    """Exhausting max_iterations without a final answer records the event."""

    def always_tool(messages: list[dict[str, Any]]) -> dict[str, Any]:
        return tool_turn()

    sink: list[dict[str, Any]] = []
    agentic_single_step(
        "CC=O", MODEL, llm_runner=always_tool, max_iterations=2, event_sink=sink
    )
    assert len(sink) == 1
    assert sink[0]["kind"] == "no_final_answer"


def test_event_sink_empty_on_success() -> None:
    """A successful final answer records no event."""
    model = ScriptedModel([final_turn()])
    sink: list[dict[str, Any]] = []
    result = agentic_single_step("CC=O", MODEL, llm_runner=model, event_sink=sink)
    assert result[0]  # non-empty pathways
    assert sink == []


def test_orchestrator_raises_not_implemented() -> None:
    """The top-level orchestrator is a scaffold and raises when invoked."""
    with pytest.raises(NotImplementedError, match="orchestrator"):
        agentic_orchestrator("CCO", MODEL)
