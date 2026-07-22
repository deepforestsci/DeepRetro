"""Tests for lossless assistant-message history preservation."""

from __future__ import annotations

import ast
import inspect
from typing import Any

import pytest

from deepretro.agents import loop, message_history
from deepretro.agents.message_history import (
    MessagePreservationError,
    append_assistant_message,
    snapshot_assistant_message,
)


class ProviderMessage:
    """Small provider-message stand-in exposing ``model_dump``."""

    def __init__(self, payload: dict[str, Any]) -> None:
        self.payload = payload

    def model_dump(self) -> dict[str, Any]:
        """Return the provider payload."""
        return self.payload


def test_snapshot_preserves_provider_fields_and_is_independent() -> None:
    """Provider-specific nested fields survive in an independent snapshot."""
    payload = {
        "role": "assistant",
        "content": None,
        "reasoning_content": {"signature": "signed-token"},
    }

    snapshot = snapshot_assistant_message(ProviderMessage(payload))
    payload["reasoning_content"]["signature"] = "changed"

    assert snapshot["reasoning_content"] == {"signature": "signed-token"}


def test_append_assistant_message_adds_and_returns_snapshot() -> None:
    """Appending stores the exact snapshot returned to the caller."""
    history: list[dict[str, Any]] = []

    snapshot = append_assistant_message(
        history,
        {"role": "assistant", "content": "done", "provider_field": 7},
    )

    assert history == [snapshot]
    assert snapshot["provider_field"] == 7


def test_snapshot_rejects_non_assistant_message() -> None:
    """Only assistant turns may cross the preservation boundary."""
    with pytest.raises(MessagePreservationError, match="Expected an assistant"):
        snapshot_assistant_message({"role": "user", "content": "hello"})


@pytest.mark.parametrize("module", [message_history, loop])
def test_every_agent_function_docstring_has_an_example(module: Any) -> None:
    """Every function in the PR modules includes a concrete usage example."""
    tree = ast.parse(inspect.getsource(module))
    missing = [
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and "Examples\n" not in (ast.get_docstring(node) or "")
    ]

    assert missing == []
