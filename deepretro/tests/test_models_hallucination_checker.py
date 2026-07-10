"""Tests for the deterministic hallucination checker and its tool wiring.

Named ``test_models_hallucination_checker`` to stay distinct from
``test_hallucination_checker`` (the heuristic algorithm's tests). These tests
run fully offline: the agent loop is driven by a scripted ``llm_runner`` so no
network is touched.
"""

from __future__ import annotations

import json
from typing import Any

import pytest

from deepretro.agents.loop import agentic_single_step
from deepretro.models.hallucination_checker import make_hallucination_checker

MODEL = "openai/gpt-4o-mini"
FINAL_ANSWER = (
    '<json>{"data": [["CCO"]], "explanation": ["reduce"], '
    '"confidence_scores": [0.8]}</json>'
)

# Sample inputs whose seeded-hash verdicts cover both branches (see the values
# computed for seed=0): a kept step and a flagged step.
KEEP_PRODUCT, KEEP_PATHWAYS = "CC=O", [["CCO"]]
FLAG_PRODUCT, FLAG_PATHWAYS = "CCO", [["CC", "O"]]


def test_checker_matches_contract_shape() -> None:
    """The checker returns a ``(status, kept)`` tuple."""
    checker = make_hallucination_checker(seed=0)
    status, kept = checker(KEEP_PRODUCT, KEEP_PATHWAYS)
    assert isinstance(status, int)
    assert isinstance(kept, list)


def test_checker_is_deterministic() -> None:
    """Same (product, pathways, seed) always yields the same verdict."""
    checker = make_hallucination_checker(seed=0)
    first = checker(FLAG_PRODUCT, FLAG_PATHWAYS)
    second = checker(FLAG_PRODUCT, FLAG_PATHWAYS)
    assert first == second


def test_seed_changes_are_reproducible() -> None:
    """A different seed is itself deterministic across calls."""
    checker_a = make_hallucination_checker(seed=7)
    checker_b = make_hallucination_checker(seed=7)
    assert checker_a(KEEP_PRODUCT, KEEP_PATHWAYS) == checker_b(
        KEEP_PRODUCT, KEEP_PATHWAYS
    )


def test_keep_verdict_returns_pathways() -> None:
    """A keep verdict maps to ``(200, pathways)``."""
    checker = make_hallucination_checker(seed=0)
    assert checker(KEEP_PRODUCT, KEEP_PATHWAYS) == (200, KEEP_PATHWAYS)


def test_flag_verdict_returns_empty() -> None:
    """A flag verdict maps to ``(400, [])``."""
    checker = make_hallucination_checker(seed=0)
    assert checker(FLAG_PRODUCT, FLAG_PATHWAYS) == (400, [])


def test_both_verdicts_are_reachable() -> None:
    """Across sample inputs, both keep and flag verdicts occur."""
    checker = make_hallucination_checker(seed=0)
    statuses = {
        checker("c1ccccc1", [["Brc1ccccc1"]])[0],
        checker("O", [["[OH2]"]])[0],
        checker(KEEP_PRODUCT, KEEP_PATHWAYS)[0],
        checker(FLAG_PRODUCT, FLAG_PATHWAYS)[0],
    }
    assert statuses == {200, 400}


def test_checker_prints_role_tagged_line(capsys: pytest.CaptureFixture[str]) -> None:
    """Every invocation prints a ``[hallucination-checker]`` line."""
    checker = make_hallucination_checker(seed=0)
    checker(KEEP_PRODUCT, KEEP_PATHWAYS)
    out = capsys.readouterr().out
    assert "[hallucination-checker] check_hallucination called:" in out
    assert "keep=True" in out


class ScriptedModel:
    """Injectable model that replays scripted turns and records inputs."""

    def __init__(self, turns: list[dict[str, Any]]) -> None:
        self.turns = list(turns)
        self.seen: list[list[dict[str, Any]]] = []

    def __call__(self, messages: list[dict[str, Any]]) -> dict[str, Any]:
        self.seen.append([dict(message) for message in messages])
        return self.turns[len(self.seen) - 1]


def _hallucination_tool_turn(call_id: str = "call_1") -> dict[str, Any]:
    """An assistant turn requesting a ``check_hallucination`` tool call."""
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


def test_agent_invokes_checker_and_feeds_result_back(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The agent calls ``check_hallucination``; its result is fed back first."""
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

    # Final answer parsed through to the pipeline shape.
    assert pathways == [["CCO"]]
    assert explanations == ["reduce"]
    assert confidence == [0.8]

    # The checker actually ran (its print line is visible).
    out = capsys.readouterr().out
    assert "[hallucination-checker] check_hallucination called:" in out

    # The tool result was appended before the model produced the final answer.
    second_call_messages = model.seen[1]
    tool_messages = [m for m in second_call_messages if m.get("role") == "tool"]
    assert tool_messages, "tool result message was not fed back to the model"
    assert tool_messages[0]["tool_call_id"] == "abc123"
    result = json.loads(tool_messages[0]["content"])
    assert result["is_hallucination"] is False
