"""Tests for offline defaults and the explicit integration opt-in."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

TESTS_DIR = Path(__file__).parent

POLICY_PROBE = """import socket
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
    with socket.create_server(("127.0.0.1", 0)) as server:
        socket.create_connection(server.getsockname()).close()
    Path("integration-ran").write_text("yes")
"""

E2E_FILES = {
    "test_az_e2e.py": (TESTS_DIR / "test_az_e2e.py").read_text(encoding="utf-8")
}


def _run_pytest(
    tmp_path: Path, test_files: dict[str, str], *options: str, **env: str
) -> subprocess.CompletedProcess[str]:
    """Run a fresh pytest in ``tmp_path`` with our conftest and the given files."""
    files = {
        "conftest.py": (TESTS_DIR / "conftest.py").read_text(encoding="utf-8"),
        "pytest.ini": "[pytest]\nmarkers =\n    integration: x\n    slow: x\n",
        **test_files,
    }
    for name, text in files.items():
        (tmp_path / name).write_text(text, encoding="utf-8")
    environment: dict[str, str] = {**os.environ, "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1"}
    environment.pop("AZ_TEST_MODELS_DIR", None)
    environment.update(env)
    return subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *options],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )


def test_unit_tests_cannot_open_sockets(tmp_path: Path) -> None:
    result = _run_pytest(
        tmp_path, {"test_probe.py": POLICY_PROBE}, "-m", "not integration"
    )

    assert "1 passed, 1 deselected" in result.stdout, result.stdout


def test_integration_tests_skip_without_run_integration_even_with_api_keys(
    tmp_path: Path,
) -> None:
    result = _run_pytest(
        tmp_path,
        {"test_probe.py": POLICY_PROBE},
        "-m",
        "integration",
        ANTHROPIC_API_KEY="not-a-real-key",
        OPENAI_API_KEY="not-a-real-key",
    )

    assert "1 skipped, 1 deselected" in result.stdout, result.stdout
    assert not (tmp_path / "integration-ran").exists()


def test_integration_tests_can_open_sockets_with_run_integration(
    tmp_path: Path,
) -> None:
    result = _run_pytest(tmp_path, {"test_probe.py": POLICY_PROBE}, "--run-integration")

    assert "2 passed" in result.stdout, result.stdout
    assert (tmp_path / "integration-ran").exists()


def test_unknown_marker_is_an_error_under_the_project_config(tmp_path: Path) -> None:
    (tmp_path / "test_typo.py").write_text(
        "import pytest\n\n@pytest.mark.integraton\ndef test_typo():\n    pass\n",
        encoding="utf-8",
    )
    project_config = TESTS_DIR.parent / "pyproject.toml"
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider"]
        + ["-c", str(project_config), str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )

    assert "'integraton' not found in `markers` configuration option" in result.stdout
    assert result.returncode == 2  # pytest.ExitCode.INTERRUPTED: collection error


def test_az_e2e_errors_when_models_dir_is_unset(tmp_path: Path) -> None:
    result = _run_pytest(tmp_path, E2E_FILES, "--run-integration")

    assert "2 errors" in result.stdout, result.stdout
    assert (
        "Set AZ_TEST_MODELS_DIR to a directory with downloaded models" in result.stdout
    )


def test_az_e2e_errors_when_config_yml_is_missing(tmp_path: Path) -> None:
    empty_dir = tmp_path / "models"
    empty_dir.mkdir()
    result = _run_pytest(
        tmp_path, E2E_FILES, "--run-integration", AZ_TEST_MODELS_DIR=str(empty_dir)
    )

    assert "2 errors" in result.stdout, result.stdout
    assert "AZ_TEST_MODELS_DIR has no config.yml" in result.stdout


def test_az_tests_error_when_az_extra_is_missing(tmp_path: Path) -> None:
    probe = (
        "import sys\nimport pytest\n\n"
        "sys.modules['aizynthfinder'] = None  # hide the az extra\n\n\n"
        "@pytest.mark.integration\n"
        "def test_needs_az(az_models_dir):\n    pass\n"
    )
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "config.yml").write_text("", encoding="utf-8")
    result = _run_pytest(
        tmp_path,
        {"test_probe.py": probe},
        "--run-integration",
        AZ_TEST_MODELS_DIR=str(models_dir),
    )

    assert "1 error" in result.stdout, result.stdout
    assert "Install the deepretro az extra to run AiZynthFinder tests" in result.stdout
