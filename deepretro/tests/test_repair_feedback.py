"""Offline regression tests for bounded feedback-driven proposal calls."""

from typing import Any

import pytest

from deepretro.algorithms.autosolve import AutoSolver
from deepretro.utils import llm, llm_interface


def test_feedback_reaches_prompt_and_caps_attempts() -> None:
    """Repair reasons reach the proposal without mutating default calls."""
    calls: list[dict[str, Any]] = []

    def runner(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
        calls.append(kwargs)
        return [["CCO"]], ["oxidation"], [0.9]

    solver = AutoSolver(llm_runner=runner, hallucination_mode="none")
    solver.run_llm("CC=O", rejection_feedback="Avoid the rejected self-loop.")
    assert calls[0]["max_attempts"] == 1
    assert "Avoid the rejected self-loop." in calls[0]["messages"][-1]["content"]
    solver.run_llm("CC=O")
    assert "max_attempts" not in calls[1]
    assert "messages" not in calls[1]


def test_feedback_is_not_silently_ignored_by_agent_mode() -> None:
    """An unsupported feedback path must fail before an external call."""
    solver = AutoSolver(solve_mode="single_step_agent")
    with pytest.raises(ValueError, match="pipeline"):
        solver.run_llm("CCO", rejection_feedback="retry")


def test_pipeline_budget_bounds_actual_transport_calls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An explicit pipeline cap disables both nested retry layers."""
    calls: list[dict[str, Any]] = []

    def unavailable(**kwargs: Any) -> None:
        calls.append(kwargs)
        raise RuntimeError("offline provider unavailable")

    monkeypatch.setattr(llm_interface, "completion", unavailable)
    assert llm.llm_pipeline("CCO", model="openai/gpt-4o-mini", max_attempts=2) == (
        [],
        [],
        [],
    )
    assert len(calls) == 2
    assert all(call["num_retries"] == 0 for call in calls)


@pytest.mark.parametrize("limit", [0, -1, True, 1.5])
def test_invalid_budgets_make_no_provider_call(
    monkeypatch: pytest.MonkeyPatch, limit: Any
) -> None:
    """Invalid limits are rejected rather than changing runtime unexpectedly."""

    def unexpected(**kwargs: Any) -> None:
        pytest.fail("Provider must not be called")

    monkeypatch.setattr(llm_interface, "completion", unexpected)
    with pytest.raises(ValueError, match="positive"):
        llm.llm_pipeline("CCO", max_attempts=limit)
    interface = llm_interface.create_llm_interface("openai/gpt-4o-mini")
    with pytest.raises(ValueError, match="positive"):
        interface.call(
            llm_interface.LLMRequest(
                "CCO", "openai/gpt-4o-mini", max_api_attempts=limit
            )
        )
