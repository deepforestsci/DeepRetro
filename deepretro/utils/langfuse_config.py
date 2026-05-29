"""Langfuse metadata configuration helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import cast

from deepretro.types import JSONValue


def get_langfuse_metadata(session_type: str = "default") -> dict[str, JSONValue]:
    """Return Langfuse metadata for a configured session type.

    Parameters
    ----------
    session_type : str, optional
        Session configuration key under ``config/langfuse_config.json``.

    Returns
    -------
    dict[str, JSONValue]
        Copy of the configured Langfuse metadata.

    Examples
    --------
    >>> metadata = get_langfuse_metadata("metadata")  # doctest: +SKIP
    >>> metadata["session_id"]  # doctest: +SKIP
    'metadata_agent'
    """
    config_path = None
    for directory in Path(__file__).resolve().parents:
        candidate = directory / "config" / "langfuse_config.json"
        if candidate.is_file():
            config_path = candidate
            break
    if config_path is None:
        raise FileNotFoundError("Could not locate config/langfuse_config.json")

    with config_path.open() as handle:
        config = json.load(handle)

    sessions = config.get("sessions", {})
    if session_type not in sessions:
        available = ", ".join(sorted(sessions))
        raise KeyError(
            f"Session type {session_type!r} not found in config. "
            f"Available sessions: {available}"
        )
    return cast(dict[str, JSONValue], sessions[session_type].copy())
