"""Tests for the single-step tool-calling agent loop."""

from __future__ import annotations

import json
from typing import Any

import pytest

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


def test_orchestrator_raises_not_implemented() -> None:
    """The top-level orchestrator is a scaffold and raises when invoked."""
    with pytest.raises(NotImplementedError, match="orchestrator"):
        agentic_orchestrator("CCO", MODEL)
