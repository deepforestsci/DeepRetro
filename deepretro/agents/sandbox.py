"""Hardened subprocess sandbox for agent-generated Python code.

The agent's ``run_python`` tool executes model-written code. That code is
*semi-trusted* (produced by the LLM's own reasoning, never interpolated from
untrusted Google-Sheet / molecule-file text), so the sandbox aims to make the
common failure modes safe rather than to be an escape-proof jail:

* **Secret scrubbing (primary guarantee).** The child runs with a minimal
  allow-list environment, so API keys in the parent process (``ANTHROPIC_API_KEY``
  and friends) can never be read or exfiltrated by generated code.
* **Resource caps.** POSIX ``RLIMIT_AS`` / ``RLIMIT_CPU`` / ``RLIMIT_NPROC`` /
  ``RLIMIT_FSIZE`` bound memory, CPU, fork-bombs, and file writes.
* **Wall-clock timeout.** The child process group is killed on expiry.
* **Best-effort network isolation.** On Linux with ``unshare`` available the
  child runs in a fresh, network-less namespace.

ponytail: subprocess + rlimits + secret-scrub is adequate for semi-trusted
code; it is NOT namespace isolation on macOS (``RLIMIT_AS`` unenforced, no
``unshare``). Upgrade path: a ``PodmanSandbox`` (``--network=none``, read-only
rootfs) implementing the same :class:`Sandbox` protocol.
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

    Attributes
    ----------
    ok : bool
        ``True`` when the code exited with status 0 within limits.
    stdout : str
        Captured standard output.
    stderr : str
        Captured standard error (traceback on failure).
    error : str
        Sandbox-level error (e.g. a timeout message); empty on success.

    Examples
    --------
    >>> SandboxResult(ok=True, stdout="hi\\n", stderr="", error="").ok
    True
    """

    ok: bool
    stdout: str
    stderr: str
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
        When ``False`` (default), the child is run under ``unshare -rn`` if
        available (Linux) to remove network access.

    Examples
    --------
    >>> SubprocessSandbox().run("print(2 + 2)").stdout.strip()
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
        >>> SubprocessSandbox().run("print('ok')").ok
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
    """Kill the child's process group so no orphaned children survive."""
    try:
        os.killpg(os.getpgid(pid), 9)
    except (ProcessLookupError, PermissionError, OSError):  # pragma: no cover
        pass
