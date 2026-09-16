"""Tests for the single-step tool-calling agent loop."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import deepretro.agents.loop as agent_loop
from deepretro.agents.loop import agentic_orchestrator, agentic_single_step
from deepretro.utils.llm_trace import LOG_FILENAME, molecule_trace

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


@pytest.mark.parametrize("backend", ["structured", "sandbox"])
def test_mask_reason_restore_workflow(backend: str) -> None:
    """Tools and guidance reach the model; only restored precursors escape."""
    calls = 0

    def model(messages: list[dict[str, Any]]) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        if calls == 1:
            assert "handle_protection" in messages[0]["content"]
            assert "handle_deprotection" in messages[0]["content"]
            return tool_turn("handle_protection", {"smiles": "COc1ccc(C=O)cc1"})
        result = json.loads(messages[-1]["content"])
        if calls == 2:
            assert result["matched"]
            assert "exposed core" in result["guidance"]
            marker = result["groups"][0]["marker"]
            return tool_turn(
                "handle_deprotection",
                {"smiles": f"{marker}Oc1ccc(CO)cc1", "mask_id": result["mask_id"]},
                call_id="restore",
            )
        return final_turn(
            "<json>"
            + json.dumps(
                {
                    "data": [[result["smiles"]]],
                    "explanation": ["Transform the core with methoxy intact"],
                    "confidence_scores": [0.8],
                }
            )
            + "</json>"
        )

    pathways, _, _ = agentic_single_step(
        "COc1ccc(C=O)cc1", MODEL, llm_runner=model, tool_backend=backend
    )
    assert pathways == [["COc1ccc(CO)cc1"]]
    assert calls == 3


def test_unresolved_masked_final_answer_requests_repair() -> None:
    bad = '<json>{"data": [["CO[*:1]"]], "explanation": ["x"], "confidence_scores": [0.8]}</json>'
    model = ScriptedModel([final_turn(bad), final_turn()])
    assert agentic_single_step("CC=O", MODEL, llm_runner=model)[0] == [["CCO"]]
    assert "unresolved dummy atoms" in model.seen[1][-1]["content"]


@pytest.mark.parametrize("backend", ["structured", "sandbox"])
def test_agent_can_propose_chemical_deprotection(backend: str) -> None:
    """The model tests a protected precursor and receives its forward product."""
    precursor = "CCNC(=O)OC(C)(C)C"
    calls = 0

    def model(messages: list[dict[str, Any]]) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        if calls == 1:
            assert "mode='propose'" in messages[0]["content"]
            return tool_turn(
                "handle_deprotection",
                {"smiles": precursor, "mode": "propose", "groups": ["Boc"]},
            )
        result = json.loads(messages[-1]["content"])
        assert result["direction"] == "forward"
        assert result["candidates"][0]["product_smiles"] == "CCN"
        return final_turn(
            "<json>"
            + json.dumps(
                {
                    "data": [[precursor]],
                    "explanation": ["Candidate Boc removal"],
                    "confidence_scores": [0.5],
                }
            )
            + "</json>"
        )

    assert agentic_single_step("CCN", MODEL, llm_runner=model, tool_backend=backend)[
        0
    ] == [[precursor]]


def test_unresolved_masks_exhaust_budget_without_returning_pathway() -> None:
    bad = '<json>{"data": [["CO[*:1]"]], "explanation": ["x"], "confidence_scores": [0.8]}</json>'
    events: list[dict[str, Any]] = []
    result = agentic_single_step(
        "COC",
        MODEL,
        llm_runner=lambda _: final_turn(bad),
        max_iterations=2,
        event_sink=events,
    )
    assert result == ([], [], [])
    assert events[0]["kind"] == "no_final_answer"


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


class TestIterationBudget:
    """Per-node agent iteration budget from carbon count and depth."""

    def test_large_molecule_at_root_hits_max_cap(self) -> None:
        assert agent_loop.iteration_budget("C" * 20, depth=0) == 15

    def test_budget_decays_with_depth(self) -> None:
        # 20 carbons: 20, 15, 11.25, 8.44, 6.33, 4.75 -> clamped to [5, 15]
        budgets = [agent_loop.iteration_budget("C" * 20, depth=d) for d in range(6)]
        assert budgets == [15, 15, 11, 8, 6, 5]

    def test_half_rounds_up_not_to_even(self) -> None:
        # 10 carbons at depth 1 -> 7.5 -> 8 (not banker's rounding to 8/7 mix)
        assert agent_loop.iteration_budget("C" * 10, depth=1) == 8
        # 6 carbons at depth 1 -> 4.5 -> 5 after rounding, also the floor
        assert agent_loop.iteration_budget("C" * 6, depth=1) == 5

    def test_small_molecule_gets_min_cap(self) -> None:
        assert agent_loop.iteration_budget("CCO", depth=0) == 5

    def test_no_carbon_gets_min_cap(self) -> None:
        assert agent_loop.iteration_budget("O=S(=O)(O)O", depth=0) == 5

    def test_unparseable_gets_min_cap(self) -> None:
        assert agent_loop.iteration_budget("not-a-smiles", depth=0) == 5

    def test_custom_bounds_and_decay(self) -> None:
        budget = agent_loop.iteration_budget(
            "C" * 10, depth=1, min_iterations=2, max_iterations=8, decay=0.5
        )
        assert budget == 5

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"min_iterations": 0},
            {"min_iterations": 6, "max_iterations": 5},
            {"decay": 0.0},
            {"decay": 1.5},
            {"depth": -1},
        ],
    )
    def test_rejects_bad_parameters(self, kwargs: dict[str, Any]) -> None:
        params: dict[str, Any] = {"depth": 0}
        params.update(kwargs)
        with pytest.raises(ValueError):
            agent_loop.iteration_budget("CCO", **params)
# ---------------------------------------------------------------------------
# Per-molecule LLM call logging (default litellm-backed model call)
# ---------------------------------------------------------------------------


class FakeCompletionMessage:
    """Assistant message stand-in returned by the patched ``litellm.completion``."""

    def __init__(
        self,
        content: str | None,
        tool_calls: list[dict[str, Any]] | None = None,
    ) -> None:
        self.content = content
        self.tool_calls = tool_calls

    def model_dump(self) -> dict[str, Any]:
        dumped: dict[str, Any] = {"role": "assistant", "content": self.content}
        if self.tool_calls is not None:
            dumped["tool_calls"] = self.tool_calls
        return dumped


class FakeCompletionResponse:
    """Minimal LiteLLM ``ModelResponse`` stand-in."""

    def __init__(
        self,
        content: str | None,
        tool_calls: list[dict[str, Any]] | None = None,
    ) -> None:
        self.choices = [FakeChoice(FakeCompletionMessage(content, tool_calls))]
        self.usage = {"prompt_tokens": 3, "completion_tokens": 5, "total_tokens": 8}


class FakeChoice:
    """LiteLLM choice stand-in exposing one assistant message."""

    def __init__(self, message: FakeCompletionMessage) -> None:
        self.message = message


def patch_completion(
    monkeypatch: pytest.MonkeyPatch, responses: list[Any]
) -> list[dict[str, Any]]:
    """Patch ``litellm.completion`` with scripted responses; return seen params."""
    import litellm

    seen: list[dict[str, Any]] = []

    def fake_completion(**params: Any) -> Any:
        seen.append(params)
        result = responses[len(seen) - 1]
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(litellm, "completion", fake_completion)
    return seen


def read_log(log_dir: Path) -> list[dict[str, Any]]:
    """Read the JSONL call log written under *log_dir*."""
    path = log_dir / LOG_FILENAME
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def test_default_model_call_records_each_iteration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Each agent turn writes one record carrying its 1-based iteration."""
    tool_answer = FakeCompletionResponse(
        None,
        [
            {
                "id": "call_1",
                "type": "function",
                "function": {
                    "name": "validate_smiles",
                    "arguments": json.dumps({"smiles": "CCO"}),
                },
            }
        ],
    )
    seen = patch_completion(
        monkeypatch, [tool_answer, FakeCompletionResponse(FINAL_ANSWER)]
    )

    with molecule_trace("CC=O", log_dir=tmp_path, session_id="sess-1"):
        result = agentic_single_step("CC=O", MODEL)

    assert result == ([["CCO"]], ["reduce"], [0.8])
    records = read_log(tmp_path)
    assert [record["iteration"] for record in records] == [1, 2]
    assert {record["stage"] for record in records} == {"retrosynthesis_agent"}
    assert records[0]["tool_calls"][0]["id"] == "call_1"
    assert records[1]["response"] == FINAL_ANSWER
    assert records[1]["usage"] == {
        "prompt_tokens": 3,
        "completion_tokens": 5,
        "total_tokens": 8,
    }
    assert seen[0]["metadata"]["session_id"] == "sess-1"
    assert seen[0]["metadata"]["generation_name"] == "retrosynthesis_agent"
    assert seen[0]["metadata"]["task"] == "retrosynthesis_agent"


def test_default_model_call_records_failures_and_reraises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A provider error is recorded and then propagated unchanged."""
    patch_completion(monkeypatch, [RuntimeError("provider down")])

    with molecule_trace("CC=O", log_dir=tmp_path, session_id="sess-1"):
        with pytest.raises(RuntimeError, match="provider down"):
            agentic_single_step("CC=O", MODEL)

    record = read_log(tmp_path)[0]
    assert record["error"] == "provider down"
    assert record["response"] is None


def test_default_model_call_without_a_trace_writes_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Outside a molecule trace no log file is created and no ids are sent."""
    seen = patch_completion(monkeypatch, [FakeCompletionResponse(FINAL_ANSWER)])

    agentic_single_step("CC=O", MODEL)

    assert not (tmp_path / LOG_FILENAME).exists()
    assert "session_id" not in seen[0]["metadata"]
    assert seen[0]["metadata"]["task"] == "retrosynthesis_agent"
