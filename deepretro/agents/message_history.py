"""Lossless in-memory preservation for provider assistant messages.

Assistant messages can contain signed, encrypted, or provider-specific
reasoning fields needed by a later tool-calling turn.  This module treats the
complete serialized message as opaque protocol data: it validates only the
outer role and returns an isolated snapshot without filtering nested fields.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from typing import Any


class MessagePreservationError(RuntimeError):
    """Raised when an assistant message cannot be preserved losslessly."""


def _serialize_message(message: Any) -> Mapping[str, Any]:
    """Return a complete mapping representation without inspecting its data."""
    if isinstance(message, Mapping):
        return message

    try:
        model_dump = getattr(message, "model_dump", None)
    except Exception:
        raise MessagePreservationError(
            "Assistant message serialization failed"
        ) from None
    if not callable(model_dump):
        raise MessagePreservationError("Unsupported assistant message object")
    try:
        serialized = model_dump()
    except Exception:
        raise MessagePreservationError(
            "Assistant message serialization failed"
        ) from None
    if not isinstance(serialized, Mapping):
        raise MessagePreservationError(
            "Assistant message serialization did not return a mapping"
        )
    return serialized


def snapshot_assistant_message(message: Any) -> dict[str, Any]:
    """Return an isolated, complete Python snapshot of an assistant message.

    Parameters
    ----------
    message : mapping or Pydantic-compatible object
        Serialized assistant mapping or an object exposing ``model_dump()``.

    Returns
    -------
    dict[str, Any]
        A recursively independent snapshot containing every serialized field.

    Raises
    ------
    MessagePreservationError
        If the object cannot be serialized, is not an assistant turn, or
        cannot be copied. Error text never includes response data.
    """
    serialized = _serialize_message(message)
    try:
        is_assistant = serialized.get("role") == "assistant"
    except Exception:
        raise MessagePreservationError(
            "Assistant message preservation failed"
        ) from None
    if not is_assistant:
        raise MessagePreservationError("Expected an assistant message")

    try:
        snapshot = deepcopy(dict(serialized))
    except Exception:
        raise MessagePreservationError(
            "Assistant message preservation failed"
        ) from None
    return snapshot


def append_assistant_message(
    history: list[dict[str, Any]],
    message: Any,
) -> dict[str, Any]:
    """Snapshot one assistant message, append it, and return the snapshot."""
    snapshot = snapshot_assistant_message(message)
    history.append(snapshot)
    return snapshot
