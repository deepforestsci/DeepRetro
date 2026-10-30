"""Configuration lookup uses checkout paths without contacting Langfuse."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from deepretro.utils import langfuse_config


@pytest.fixture
def config_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Place a fake module below a temporary checkout configuration."""
    monkeypatch.setattr(
        langfuse_config, "__file__", str(tmp_path / "deepretro" / "utils" / "config.py")
    )
    directory = tmp_path / "config"
    directory.mkdir()
    return directory / "langfuse_config.json"


def test_config_reads_named_session(config_path: Path) -> None:
    """Session metadata is read from the nearest ancestor configuration."""
    config_path.write_text(
        json.dumps({"sessions": {"metadata": {"session_id": "metadata_agent"}}}),
        encoding="utf-8",
    )
    metadata = langfuse_config.get_langfuse_metadata("metadata")
    assert metadata == {"session_id": "metadata_agent"}
    metadata["session_id"] = "changed"
    assert (
        langfuse_config.get_langfuse_metadata("metadata")["session_id"]
        == "metadata_agent"
    )


def test_config_reports_missing_session(config_path: Path) -> None:
    """The error identifies both the missing key and available session names."""
    config_path.write_text('{"sessions": {"default": {}}}', encoding="utf-8")
    with pytest.raises(KeyError, match="'missing'.*Available sessions: default"):
        langfuse_config.get_langfuse_metadata("missing")


def test_config_reports_missing_file(config_path: Path) -> None:
    """A source checkout without a configuration has an actionable error."""
    with pytest.raises(FileNotFoundError, match="config/langfuse_config.json"):
        langfuse_config.get_langfuse_metadata()


def test_config_reports_invalid_json(config_path: Path) -> None:
    """Malformed configuration is not silently replaced with defaults."""
    config_path.write_text("{broken", encoding="utf-8")
    with pytest.raises(json.JSONDecodeError):
        langfuse_config.get_langfuse_metadata()
