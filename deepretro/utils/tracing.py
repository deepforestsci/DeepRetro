"""Optional Langfuse tracing shim.

Langfuse is a *degradable* dependency: the package declares it, but a given
runtime may not have it installed or may run without ``LANGFUSE_*`` env vars
configured. This module exposes an :func:`observe` decorator and a
:func:`flush_traces` helper that use Langfuse when it is importable and are safe
no-ops otherwise, so callers never have to branch on whether tracing is present.

Examples
--------
>>> @observe(name="unit")
... def add(a: int, b: int) -> int:
...     return a + b
>>> add(2, 3)
5
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar

import structlog

logger = structlog.get_logger(__name__)

F = TypeVar("F", bound=Callable[..., Any])

try:  # langfuse<3 exposes the decorator API under ``langfuse.decorators``
    from langfuse.decorators import langfuse_context, observe

    _TRACING_AVAILABLE = True
except Exception:  # pragma: no cover - exercised only when langfuse is absent
    langfuse_context = None

    def observe(*decorator_args: Any, **decorator_kwargs: Any) -> Any:
        """No-op fallback for langfuse ``observe`` when it cannot be imported.

        Supports both ``@observe`` and ``@observe(name=...)`` usage and returns
        the wrapped function unchanged.
        """
        if decorator_args and callable(decorator_args[0]):
            return decorator_args[0]

        def _decorate(func: F) -> F:
            return func

        return _decorate

    _TRACING_AVAILABLE = False


def tracing_available() -> bool:
    """Return whether Langfuse tracing is importable in this environment.

    Returns
    -------
    bool
        ``True`` when the langfuse decorator API imported successfully.

    Examples
    --------
    >>> isinstance(tracing_available(), bool)
    True
    """
    return _TRACING_AVAILABLE


def flush_traces() -> None:
    """Flush pending Langfuse traces when tracing is available; else do nothing.

    Examples
    --------
    >>> flush_traces()  # safe whether or not langfuse is installed
    """
    if not _TRACING_AVAILABLE or langfuse_context is None:
        return
    try:
        langfuse_context.flush()
    except Exception as exc:  # pragma: no cover - network/config dependent
        logger.info("tracing.flush_unavailable", error=str(exc))
