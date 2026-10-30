"""Run agent-generated Python in a temporary subprocess.

The child receives a small allow-list of environment variables. Provider API
keys are not inherited. POSIX resource limits and a wall-clock timeout reduce
accidental resource use. On Linux, an available ``unshare`` command requests a
separate network namespace.

This is not a security boundary for untrusted code. The child keeps the user's
filesystem permissions, including access to files under HOME. On macOS there
is no network namespace and some resource limits are not enforced. Supply a
separately isolated backend implementing :class:`Sandbox` when stronger
isolation is required.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import structlog

logger = structlog.get_logger(__name__)

# Environment variables the child is allowed to inherit. Everything else
# (notably provider API keys) is dropped.
_ENV_ALLOWLIST = ("PATH", "LANG", "LC_ALL", "LC_CTYPE", "HOME", "TMPDIR", "SYSTEMROOT")

# Floor for the address-space limit so RDKit/deepretro imports do not
# MemoryError on start-up regardless of a small requested ``mem_mb``.
_MIN_ADDRESS_SPACE_MB = 1024


@dataclass(frozen=True)
class SandboxResult:
    """Outcome of one sandboxed execution.

    Examples
    --------
    >>> SandboxResult(ok=True, stdout="hi\\n", stderr="", error="").ok
    True
    """

    #: ``True`` when the code exited with status 0 within limits.
    ok: bool
    #: Captured standard output.
    stdout: str
    #: Captured standard error (traceback on failure).
    stderr: str
    #: Sandbox-level error (e.g. a timeout message); empty on success.
    error: str


@runtime_checkable
class Sandbox(Protocol):
    """Protocol for code-execution backends used by the ``run_python`` tool."""

    def run(
        self, code: str, *, timeout_s: float = ..., mem_mb: int = ...
    ) -> SandboxResult:
        """Execute *code* and return a :class:`SandboxResult`."""
        ...


def _build_child_env() -> dict[str, str]:
    """Return a minimal environment with provider secrets removed.

    Returns
    -------
    dict[str, str]
        Only the allow-listed variables from the parent environment.

    Examples
    --------
    >>> "ANTHROPIC_API_KEY" in _build_child_env()
    False
    """
    env = {key: os.environ[key] for key in _ENV_ALLOWLIST if key in os.environ}
    # Force unbuffered output so captured stdout is complete even on kill.
    env["PYTHONUNBUFFERED"] = "1"
    return env


def _apply_limits(mem_mb: int, cpu_seconds: int) -> None:
    """Start a new session and apply POSIX resource limits (best effort).

    Runs in the child via ``preexec_fn``. Each limit is guarded because some
    are unsupported or unenforced on certain platforms (e.g. macOS).
    """
    import resource

    os.setsid()
    mem_bytes = max(mem_mb, _MIN_ADDRESS_SPACE_MB) * 1024 * 1024
    limits = [
        (resource.RLIMIT_AS, mem_bytes),
        (resource.RLIMIT_CPU, cpu_seconds),
        (resource.RLIMIT_FSIZE, 64 * 1024 * 1024),
    ]
    if hasattr(resource, "RLIMIT_NPROC"):
        limits.append((resource.RLIMIT_NPROC, 64))
    for res, limit in limits:
        try:
            resource.setrlimit(res, (limit, limit))
        except (ValueError, OSError):  # pragma: no cover - platform dependent
            pass


class SubprocessSandbox:
    """Run Python code in a resource-limited, secret-scrubbed subprocess.

    Parameters
    ----------
    python_executable : str, optional
        Interpreter used for the child. Defaults to the current interpreter,
        so installed packages (RDKit, deepretro) remain importable.
    allow_network : bool, optional
        When ``False`` (default), use ``unshare -rn`` on Linux if installed.
        Other platforms still allow network access. A failing ``unshare``
        command is reported as a failed execution; it is not retried without
        isolation.

    Examples
    --------
    This trusted example bypasses the optional Linux namespace. It does not
    make network requests.

    >>> SubprocessSandbox(allow_network=True).run("print(2 + 2)").stdout.strip()
    '4'
    """

    def __init__(
        self,
        python_executable: str | None = None,
        *,
        allow_network: bool = False,
    ) -> None:
        """Store the interpreter path and network policy."""
        self.python_executable = python_executable or sys.executable
        self.allow_network = allow_network

    def _command(self, script_path: str) -> list[str]:
        """Build the argv, prefixing an unshare wrapper when isolating network."""
        base = [self.python_executable, "-I", script_path]
        if self.allow_network:
            return base
        unshare = shutil.which("unshare")
        if unshare is not None and sys.platform.startswith("linux"):
            return [unshare, "-rn", *base]
        return base

    def run(
        self,
        code: str,
        *,
        timeout_s: float = 10.0,
        mem_mb: int = 2048,
    ) -> SandboxResult:
        """Execute *code* and capture its output.

        Parameters
        ----------
        code : str
            Python source to execute.
        timeout_s : float, optional
            Wall-clock timeout in seconds.
        mem_mb : int, optional
            Address-space limit in MiB (floored at 1024).

        Returns
        -------
        SandboxResult
            Execution outcome.

        Examples
        --------
        >>> SubprocessSandbox(allow_network=True).run("print('ok')").ok
        True
        """
        cpu_seconds = max(1, int(timeout_s) + 1)
        preexec = (
            (lambda: _apply_limits(mem_mb, cpu_seconds)) if os.name == "posix" else None
        )
        with tempfile.TemporaryDirectory(prefix="deepretro_sandbox_") as workdir:
            script_path = os.path.join(workdir, "snippet.py")
            with open(script_path, "w", encoding="utf-8") as handle:
                handle.write(code)

            try:
                proc = subprocess.Popen(
                    self._command(script_path),
                    cwd=workdir,
                    env=_build_child_env(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    preexec_fn=preexec,
                )
            except OSError as exc:  # pragma: no cover - spawn failure
                logger.error("Sandbox failed to start", error=str(exc))
                return SandboxResult(False, "", "", f"failed to start sandbox: {exc}")

            try:
                stdout, stderr = proc.communicate(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                _kill_process_group(proc.pid)
                stdout, stderr = proc.communicate()
                return SandboxResult(
                    False,
                    stdout or "",
                    stderr or "",
                    f"execution timed out after {timeout_s}s",
                )

            ok = proc.returncode == 0
            error = "" if ok else f"exited with status {proc.returncode}"
            return SandboxResult(ok, stdout or "", stderr or "", error)


def _kill_process_group(pid: int) -> None:
    """Kill processes still in the child's process group on timeout."""
    try:
        os.killpg(os.getpgid(pid), 9)
    except (ProcessLookupError, PermissionError, OSError):  # pragma: no cover
        pass
