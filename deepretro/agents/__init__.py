"""Agentic tool-calling layer for DeepRetro retrosynthesis.

Exposes the code-execution sandbox, the tool registry, and the single-step
agent loop used by :class:`deepretro.algorithms.autosolve.AutoSolver` when
``solve_mode="single_step_agent"``.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "SandboxResult",
    "SubprocessSandbox",
    "build_tool_registry",
    "agentic_single_step",
    "agentic_orchestrator",
]


def __getattr__(name: str) -> Any:
    """Lazily resolve agent exports (PEP 562) to keep imports light."""
    if name in ("SandboxResult", "SubprocessSandbox"):
        from deepretro.agents import sandbox

        return getattr(sandbox, name)
    if name == "build_tool_registry":
        from deepretro.agents.tools import build_tool_registry

        return build_tool_registry
    if name in ("agentic_single_step", "agentic_orchestrator"):
        from deepretro.agents import loop

        return getattr(loop, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
