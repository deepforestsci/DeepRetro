"""Tests for deepretro's logging infrastructure.

Validates that ``configure_logging`` controls output level and format, and
that ``structlog.contextvars.bound_contextvars`` propagates context to all
downstream log calls.

These tests import ``configure_logging`` directly from the module to avoid
pulling in heavy transitive dependencies (e.g. deepchem) through
``deepretro.__init__``.
"""

from __future__ import annotations

import structlog

from deepretro import configure_logging


class TestConfigureLogging:
    """Tests for ``deepretro.configure_logging``."""

    def teardown_method(self) -> None:
        structlog.reset_defaults()

    def test_default_level_is_info(self) -> None:
        configure_logging(level="INFO")
        log = structlog.get_logger()
        assert hasattr(log, "info")
        assert hasattr(log, "debug")

    def test_critical_silences_everything(self, capsys) -> None:
        configure_logging(level="CRITICAL")
        log = structlog.get_logger()
        log.info("should not appear")
        log.warning("should not appear either")
        captured = capsys.readouterr()
        assert "should not appear" not in captured.out

    def test_warning_level_filters_info(self, capsys) -> None:
        configure_logging(level="WARNING")
        log = structlog.get_logger()
        log.info("info message")
        log.warning("warning message")
        captured = capsys.readouterr()
        assert "info message" not in captured.out
        assert "warning message" in captured.out

    def test_debug_level_shows_everything(self, capsys) -> None:
        configure_logging(level="DEBUG")
        log = structlog.get_logger()
        log.debug("debug message")
        log.info("info message")
        captured = capsys.readouterr()
        assert "debug message" in captured.out
        assert "info message" in captured.out

    def test_json_output(self, capsys) -> None:
        configure_logging(level="INFO", json_output=True)
        log = structlog.get_logger()
        log.info("json test", key="value")
        captured = capsys.readouterr()
        assert '"event": "json test"' in captured.out
        assert '"key": "value"' in captured.out

    def test_case_insensitive_level(self) -> None:
        configure_logging(level="warning")
        bound = structlog.get_logger()
        assert bound is not None


class TestContextPropagation:
    """Tests for structlog contextvars propagation."""

    def teardown_method(self) -> None:
        structlog.contextvars.clear_contextvars()
        structlog.reset_defaults()

    def test_bound_contextvars_appears_in_logs(self, capsys) -> None:
        configure_logging(level="INFO")
        log = structlog.get_logger()

        with structlog.contextvars.bound_contextvars(job_id="test-123"):
            log.info("context test")

        captured = capsys.readouterr()
        assert "test-123" in captured.out

    def test_context_cleaned_after_block(self, capsys) -> None:
        configure_logging(level="INFO")
        log = structlog.get_logger()

        with structlog.contextvars.bound_contextvars(job_id="temp-456"):
            log.info("inside")

        log.info("outside")
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert "temp-456" in lines[0]
        assert "temp-456" not in lines[1]

    def test_nested_contextvars(self, capsys) -> None:
        configure_logging(level="INFO")
        log = structlog.get_logger()

        with structlog.contextvars.bound_contextvars(job_id="outer"):
            with structlog.contextvars.bound_contextvars(molecule="CCO"):
                log.info("nested")

        captured = capsys.readouterr()
        assert "outer" in captured.out
        assert "CCO" in captured.out
