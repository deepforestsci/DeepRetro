"""Regression tests for V2 decisions at tool, route, and CLI boundaries."""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from deepretro import batch
from deepretro.agents.tools import build_tool_registry
from deepretro.algorithms.autosolve import AutoSolver
from deepretro.algorithms.hallucination_weights import HallucinationWeights
from deepretro.models.hallucination_checker import HallucinationChecker


def test_retained_flagged_candidate_is_not_reported_clean() -> None:
    checker = HallucinationChecker(checker_type="heuristic")
    product, reactants = "CCCCCCCCCCCCCCCC", ["C"]
    status, retained = checker(product, [reactants])
    assert status == 200
    assert retained == [reactants]
    result = build_tool_registry(hallucination_checker=checker).execute(
        "check_hallucination", {"product": product, "reactants": reactants}
    )
    assert result["is_hallucination"] is True


def test_tool_and_route_use_configured_warning_threshold() -> None:
    weights = HallucinationWeights(reject_below=100)
    checker = HallucinationChecker(checker_type="heuristic", weights=weights)
    solver = AutoSolver(hallucination_weights=weights)
    product, reactants = "CCCC", ["CC"]
    result = build_tool_registry(hallucination_checker=checker).execute(
        "check_hallucination", {"product": product, "reactants": reactants}
    )
    verdict = solver._verdict_for(reactants, product)
    assert verdict is not None
    assert verdict["flagged"] is True
    assert result["is_hallucination"] is verdict["flagged"]


def test_disabled_checker_has_no_tool_or_route_verdict() -> None:
    solver = AutoSolver(hallucination_mode="none")
    assert solver.hallucination_checker is None
    assert solver._verdict_for(["C"], "CCCC") is None
    registry = build_tool_registry(hallucination_checker=solver.hallucination_checker)
    assert "check_hallucination" not in registry.names


@pytest.mark.parametrize("mode", ["heuristic", "none"])
def test_batch_without_training_needs_no_sheet(
    monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    download = Mock(side_effect=AssertionError("unexpected download"))
    train = Mock(side_effect=AssertionError("unexpected training"))
    run = Mock()
    monkeypatch.setattr(batch, "download_sheet_csv", download)
    monkeypatch.setattr(batch, "train_hallucination_checker", train)
    monkeypatch.setattr(batch, "read_molecules", Mock(return_value=[]))
    monkeypatch.setattr(batch, "run_batch", run)
    arguments = ["--molecules", "unused.txt", "--hallucination-mode", mode]
    if mode == "none":
        arguments.extend(["--hallucination-weights", "unused-weights.json"])
    batch.main(arguments)
    download.assert_not_called()
    train.assert_not_called()
    run.assert_called_once()


def test_batch_training_requires_sheet_before_network_access(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    download = Mock(side_effect=AssertionError("unexpected download"))
    monkeypatch.setattr(batch, "download_sheet_csv", download)
    with pytest.raises(SystemExit, match="2"):
        batch.main(["--molecules", "unused.txt", "--hallucination-mode", "ml"])
    download.assert_not_called()
