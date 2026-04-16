"""Port of src.utils.langfuse_config — resolves repo config/ from any package depth."""

import json
from pathlib import Path
from typing import Any, Dict


def _langfuse_config_path() -> Path:
    here = Path(__file__).resolve()
    for d in [here, *here.parents]:
        candidate = d / "config" / "langfuse_config.json"
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Langfuse config not found: searched upward from "
        f"{here} for config/langfuse_config.json"
    )


def get_langfuse_metadata(session_type: str = "default") -> Dict[str, Any]:
    """Get Langfuse metadata from config file.

    Args:
        session_type: Type of session ("default", "metadata", "retrosynthesis")

    Returns:
        Metadata dictionary for Langfuse

    Raises:
        FileNotFoundError: If config file doesn't exist
        KeyError: If session_type doesn't exist in config
        json.JSONDecodeError: If config file is invalid JSON
    """
    config_path = _langfuse_config_path()

    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    sessions = config.get("sessions", {})

    if session_type in sessions:
        return sessions[session_type].copy()
    raise KeyError(
        f"Session type '{session_type}' not found in config. "
        f"Available sessions: {list(sessions.keys())}"
    )
