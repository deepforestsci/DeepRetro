"""Regression tests for offline defaults and explicit integration opt-in."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize("enable_integration", [False, True])
def test_integration_requires_explicit_option(
    tmp_path: Path, enable_integration: bool
) -> None:
    """Credentials alone must not activate integration tests or unit networking."""
    conftest = Path(__file__).with_name("conftest.py")
    (tmp_path / "conftest.py").write_text(conftest.read_text(), encoding="utf-8")
    (tmp_path / "pytest.ini").write_text(
        "[pytest]\nmarkers =\n    integration: external dependency\n",
        encoding="utf-8",
    )
    (tmp_path / "test_policy_probe.py").write_text(
        """import socket
from pathlib import Path
import pytest


def test_unit_network_is_blocked():
    with pytest.raises(RuntimeError, match="Network access is disabled"):
        socket.create_connection(("127.0.0.1", 1))
    with socket.socket() as connection:
        with pytest.raises(RuntimeError, match="Network access is disabled"):
            connection.connect(("127.0.0.1", 1))
        with pytest.raises(RuntimeError, match="Network access is disabled"):
            connection.connect_ex(("127.0.0.1", 1))


@pytest.mark.integration
def test_explicit_integration():
    Path("integration-ran").write_text("yes")
""",
        encoding="utf-8",
    )
    environment = {
        **os.environ,
        "ANTHROPIC_API_KEY": "not-a-real-key",
        "OPENAI_API_KEY": "not-a-real-key",
        "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
    }
    command = [sys.executable, "-m", "pytest", "-q", "--strict-markers"]
    if enable_integration:
        command.append("--run-integration")
    result = subprocess.run(
        command,
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "integration-ran").exists() is enable_integration
    expected = "2 passed" if enable_integration else "1 passed, 1 skipped"
    assert expected in result.stdout


@pytest.mark.parametrize("configured", [False, True])
def test_az_e2e_fails_when_models_are_missing(tmp_path: Path, configured: bool) -> None:
    """An explicitly requested E2E run cannot succeed by skipping missing data."""
    for name in ("conftest.py", "test_az_e2e.py"):
        source = Path(__file__).with_name(name)
        (tmp_path / name).write_text(source.read_text(), encoding="utf-8")
    (tmp_path / "pytest.ini").write_text(
        "[pytest]\nmarkers =\n    integration: external dependency\n    slow: model search\n",
        encoding="utf-8",
    )
    environment = {**os.environ, "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1"}
    environment.pop("AZ_TEST_MODELS_DIR", None)
    if configured:
        environment["AZ_TEST_MODELS_DIR"] = str(tmp_path / "missing-models")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--run-integration"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 1, result.stdout + result.stderr
    assert "1 failed" in result.stdout
    assert (
        "Missing AiZynthFinder configuration" in result.stdout
        if configured
        else "AZ E2E requires" in result.stdout
    )
