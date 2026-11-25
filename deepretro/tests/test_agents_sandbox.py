"""Tests for the subprocess code-execution sandbox."""

from __future__ import annotations

import os
import shutil
import sys

import pytest

from deepretro.agents.sandbox import SandboxResult, SubprocessSandbox


def test_runs_code_and_captures_stdout() -> None:
    """A successful run returns ok with the captured stdout."""
    result = SubprocessSandbox(allow_network=True).run("print('hello world')")
    assert isinstance(result, SandboxResult)
    assert result.ok is True
    assert "hello world" in result.stdout


def test_reports_error_on_exception() -> None:
    """An exception in the code is reported as not-ok with a traceback."""
    result = SubprocessSandbox(allow_network=True).run("raise ValueError('boom')")
    assert result.ok is False
    assert "ValueError" in result.stderr


def test_installed_packages_are_importable() -> None:
    """RDKit (an installed dependency) imports despite the scrubbed env."""
    result = SubprocessSandbox(allow_network=True).run(
        "from rdkit import Chem; print(Chem.CanonSmiles('C(O)C'))"
    )
    assert result.ok is True
    assert "CCO" in result.stdout


def test_scrubs_secrets_from_child_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """API keys in the parent environment are not visible to child code."""
    monkeypatch.setenv("ANTHROPIC_API_KEY", "super-secret-value")
    result = SubprocessSandbox(allow_network=True).run(
        "import os; print(os.environ.get('ANTHROPIC_API_KEY', 'ABSENT'))"
    )
    assert result.ok is True
    assert "super-secret-value" not in result.stdout
    assert "ABSENT" in result.stdout


def test_times_out_on_long_running_code() -> None:
    """Wall-clock timeout kills the child and reports a timeout error."""
    result = SubprocessSandbox(allow_network=True).run(
        "import time; time.sleep(30)", timeout_s=1.0
    )
    assert result.ok is False
    assert "tim" in result.error.lower()


def test_working_directory_is_ephemeral() -> None:
    """The sandbox runs in a temp cwd that is removed after the run."""
    result = SubprocessSandbox(allow_network=True).run("import os; print(os.getcwd())")
    assert result.ok is True
    workdir = result.stdout.strip().splitlines()[-1]
    assert not os.path.exists(workdir)


@pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="RLIMIT_AS is not enforced on macOS",
)
def test_memory_limit_blocks_large_allocation() -> None:
    """A large allocation is blocked by the address-space limit (Linux only)."""
    result = SubprocessSandbox(allow_network=True).run(
        "x = bytearray(4 * 1024 * 1024 * 1024)", mem_mb=256
    )
    assert result.ok is False


def test_linux_default_requests_network_isolation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The default Linux command retains the network namespace request."""
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(shutil, "which", lambda command: "/usr/bin/unshare")
    command = SubprocessSandbox()._command("example.py")
    assert command == ["/usr/bin/unshare", "-rn", sys.executable, "-I", "example.py"]


def test_failed_namespace_does_not_retry_without_isolation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unavailable namespace must fail instead of executing the payload."""
    false_command = shutil.which("false")
    if false_command is None:
        pytest.skip("Requires the POSIX false utility")
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(shutil, "which", lambda command: false_command)
    result = SubprocessSandbox().run("print('payload ran')")
    assert not result.ok
    assert "payload ran" not in result.stdout
    assert "exited with status" in result.error
