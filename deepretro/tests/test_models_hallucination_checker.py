"""Tests for the temporary checker skeleton and tool wiring."""

from __future__ import annotations

import json
from typing import Any

import pytest

from deepretro.agents.loop import agentic_single_step
from deepretro.models.hallucination_checker import make_hallucination_checker
from deepretro.models.hallucination_helpers import resolve_hallucination

MODEL = "openai/gpt-4o-mini"
FINAL_ANSWER = (
    '<json>{"data": [["CCO"]], "explanation": ["reduce"], '
    '"confidence_scores": [0.8]}</json>'
)

KEEP_PRODUCT, KEEP_PATHWAYS = "CC=O", [["CCO"]]
FLAG_PRODUCT, FLAG_PATHWAYS = "CCO", [["CC", "O"]]


def test_default_mode_is_not_registered() -> None:
    """The temporary checker is not a public filtering mode."""
    with pytest.raises(ValueError, match="hallucination_mode must be"):
        resolve_hallucination("default", None)


def test_checker_matches_contract_shape() -> None:
    checker = make_hallucination_checker(seed=0)
    status, kept = checker(KEEP_PRODUCT, KEEP_PATHWAYS)
    assert isinstance(status, int)
    assert isinstance(kept, list)


def test_checker_is_deterministic() -> None:
    checker = make_hallucination_checker(seed=0)
    first = checker(FLAG_PRODUCT, FLAG_PATHWAYS)
    second = checker(FLAG_PRODUCT, FLAG_PATHWAYS)
    assert first == second


def test_seed_changes_are_reproducible() -> None:
    checker_a = make_hallucination_checker(seed=7)
    checker_b = make_hallucination_checker(seed=7)
    assert checker_a(KEEP_PRODUCT, KEEP_PATHWAYS) == checker_b(
        KEEP_PRODUCT, KEEP_PATHWAYS
    )


def test_keep_verdict_returns_pathways() -> None:
    checker = make_hallucination_checker(seed=0)
    assert checker(KEEP_PRODUCT, KEEP_PATHWAYS) == (200, KEEP_PATHWAYS)


def test_flag_verdict_returns_empty() -> None:
    checker = make_hallucination_checker(seed=0)
    assert checker(FLAG_PRODUCT, FLAG_PATHWAYS) == (400, [])


def test_both_verdicts_are_reachable() -> None:
    checker = make_hallucination_checker(seed=0)
    statuses = {
        checker("c1ccccc1", [["Brc1ccccc1"]])[0],
        checker("O", [["[OH2]"]])[0],
        checker(KEEP_PRODUCT, KEEP_PATHWAYS)[0],
        checker(FLAG_PRODUCT, FLAG_PATHWAYS)[0],
    }
    assert statuses == {200, 400}


class ScriptedModel:
    def __init__(self, turns: list[dict[str, Any]]) -> None:
        self.turns = list(turns)
        self.seen: list[list[dict[str, Any]]] = []

    def __call__(self, messages: list[dict[str, Any]]) -> dict[str, Any]:
        self.seen.append([dict(message) for message in messages])
        return self.turns[len(self.seen) - 1]


def _hallucination_tool_turn(call_id: str = "call_1") -> dict[str, Any]:
    return {
        "role": "assistant",
        "content": None,
        "tool_calls": [
            {
                "id": call_id,
                "type": "function",
                "function": {
                    "name": "check_hallucination",
                    "arguments": json.dumps(
                        {"product": KEEP_PRODUCT, "reactants": ["CCO"]}
                    ),
                },
            }
        ],
    }


def test_agent_invokes_checker_and_feeds_result_back() -> None:
    checker = make_hallucination_checker(seed=0)
    model = ScriptedModel(
        [
            _hallucination_tool_turn(call_id="abc123"),
            {"role": "assistant", "content": FINAL_ANSWER},
        ]
    )

    pathways, explanations, confidence = agentic_single_step(
        KEEP_PRODUCT,
        MODEL,
        hallucination_checker=checker,
        llm_runner=model,
        max_iterations=6,
    )

    assert pathways == [["CCO"]]
    assert explanations == ["reduce"]
    assert confidence == [0.8]

    second_call_messages = model.seen[1]
    tool_messages = [m for m in second_call_messages if m.get("role") == "tool"]
    assert tool_messages, "tool result message was not fed back to the model"
    assert tool_messages[0]["tool_call_id"] == "abc123"
    result = json.loads(tool_messages[0]["content"])
    assert result["is_hallucination"] is False
