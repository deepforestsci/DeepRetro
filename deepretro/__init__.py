"""
deepretro — retrosynthesis ML utilities.

Provides DeepChem-compatible featurizers for reaction-step data.
"""

from __future__ import annotations

import logging

import structlog

from deepretro.featurizers.reactionstep import ReactionStepFeaturizer

__all__ = ["ReactionStepFeaturizer", "configure_logging"]


def configure_logging(
    level: str = "INFO",
    json_output: bool = False,
) -> None:
    """Configure deepretro's logging output.

    This is an **optional** convenience function.  If never called, structlog's
    defaults apply (all levels shown, colored console output).  Applications
    that already call ``structlog.configure()`` elsewhere do not need this.

    Parameters
    ----------
    level : str
        Log level name (``"DEBUG"``, ``"INFO"``, ``"WARNING"``, ``"ERROR"``,
        ``"CRITICAL"``).  Use ``"CRITICAL"`` to effectively silence all
        logging.
    json_output : bool
        If ``True``, render logs as JSON (useful for production / log
        aggregation).  If ``False``, use colored console output.
    """
    renderer: structlog.types.Processor = (structlog.processors.JSONRenderer()
                                           if json_output else
                                           structlog.dev.ConsoleRenderer())
    structlog.configure(
        processors=[
            structlog.contextvars.merge_contextvars,
            structlog.processors.add_log_level,
            structlog.processors.StackInfoRenderer(),
            structlog.dev.set_exc_info,
            structlog.processors.TimeStamper(fmt="%Y-%m-%d %H:%M:%S",
                                             utc=False),
            renderer,
        ],
        wrapper_class=structlog.make_filtering_bound_logger(
            getattr(logging, level.upper(), logging.INFO)),
        context_class=dict,
        logger_factory=structlog.PrintLoggerFactory(),
        cache_logger_on_first_use=False,
    )
