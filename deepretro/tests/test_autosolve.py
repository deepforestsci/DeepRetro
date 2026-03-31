"""Unit tests for deepretro.utils.autosolve.

Tests the autosolve and autosolve_async wrappers using mocked pipeline
calls so the tests run without AiZynthFinder models, LLM API keys, or
any other heavy infrastructure.
"""

from __future__ import annotations

import asyncio
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from deepretro.utils.autosolve import (
    _build_ml_checker,
    _resolve_hallucination_args,
    autosolve,
    autosolve_async,
)

BENZENE = "c1ccccc1"
ETHANOL = "CCO"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

FAKE_RESULT = {
    "dependencies": {"1": []},
    "steps": [
        {
            "step": "1",
            "reactants": [{"smiles": ETHANOL}],
            "reagents": [],
            "products": [{"smiles": BENZENE}],
            "conditions": [],
            "reactionmetrics": [{"confidenceestimate": 0.85}],
        }
    ],
}


@pytest.fixture()
def mock_pipeline():
    """Patch ``_get_pipeline`` so tests never touch src/ or aizynthfinder."""
    mock_run = MagicMock(return_value=FAKE_RESULT)
    with patch(
        "deepretro.utils.autosolve._get_pipeline",
        return_value=mock_run,
    ):
        yield mock_run


class TestAutosolve:
    """Tests for ``autosolve``."""

    def test_returns_result_dict(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        result = autosolve(BENZENE)
        assert result == FAKE_RESULT
        mock_run.assert_called_once()

    def test_passes_all_params_to_pipeline(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        autosolve(
            ASPIRIN,
            llm="test-model",
            az_model="USPTO",
            stability_flag="False",
            hallucination_mode="none",
            use_protecting_group_feature=True,
        )
        mock_run.assert_called_once_with(
            ASPIRIN,
            llm="test-model",
            az_model="USPTO",
            stability_flag="False",
            hallucination_check="False",
            use_protecting_group_feature=True,
            hallucination_checker_fn=None,
        )

    def test_does_not_write_to_disk(self, mock_pipeline, tmp_path: Path) -> None:
        result = autosolve(BENZENE)
        assert result == FAKE_RESULT
        assert list(tmp_path.iterdir()) == []


class TestAutosolveAsync:
    """Tests for ``autosolve_async``."""

    def test_returns_one_result_per_input(self, mock_pipeline) -> None:
        results = asyncio.run(autosolve_async([BENZENE, ETHANOL]))
        assert len(results) == 2
        assert all(r == FAKE_RESULT for r in results)

    def test_empty_list_returns_empty(self, mock_pipeline) -> None:
        results = asyncio.run(autosolve_async([]))
        assert results == []

    def test_failed_molecule_returns_error_dict(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        mock_run.side_effect = RuntimeError("LLM down")
        results = asyncio.run(autosolve_async([BENZENE]))
        assert len(results) == 1
        assert "error" in results[0]
        assert results[0]["smiles"] == BENZENE
        assert "LLM down" in results[0]["error"]

    def test_partial_failure_does_not_crash_batch(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        mock_run.side_effect = [FAKE_RESULT, RuntimeError("fail"), FAKE_RESULT]
        results = asyncio.run(autosolve_async([BENZENE, ETHANOL, ASPIRIN]))
        assert len(results) == 3
        assert results[0] == FAKE_RESULT
        assert "error" in results[1]
        assert results[2] == FAKE_RESULT

    def test_does_not_write_to_disk(self, mock_pipeline, tmp_path: Path) -> None:
        asyncio.run(autosolve_async([BENZENE, ETHANOL]))
        assert list(tmp_path.iterdir()) == []

    def test_respects_max_concurrent(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        molecules = [f"C{'C' * i}O" for i in range(10)]
        results = asyncio.run(autosolve_async(molecules, max_concurrent=2))
        assert len(results) == 10
        assert mock_run.call_count == 10


class TestResolveHallucinationArgs:
    """Tests for ``_resolve_hallucination_args``."""

    def test_heuristic_mode(self) -> None:
        h_check, h_fn = _resolve_hallucination_args("heuristic", None)
        assert h_check == "True"
        assert h_fn is None

    def test_none_mode(self) -> None:
        h_check, h_fn = _resolve_hallucination_args("none", None)
        assert h_check == "False"
        assert h_fn is None

    def test_invalid_mode_raises(self) -> None:
        with pytest.raises(ValueError, match="hallucination_mode must be"):
            _resolve_hallucination_args("invalid", None)

    def test_ml_mode_without_classifier_raises(self) -> None:
        with pytest.raises(ValueError, match="requires hallucination_classifier"):
            _resolve_hallucination_args("ml", None)

    def test_ml_mode_with_classifier_returns_callable(self) -> None:
        mock_clf = MagicMock()
        with patch(
            "deepretro.utils.autosolve._load_classifier", return_value=mock_clf
        ), patch(
            "deepretro.utils.autosolve._build_ml_checker",
            return_value=lambda p, pw: (200, pw),
        ) as mock_build:
            h_check, h_fn = _resolve_hallucination_args("ml", mock_clf)

        assert h_check == "True"
        assert h_fn is not None
        mock_build.assert_called_once_with(mock_clf)


class TestBuildMlChecker:
    """Tests for ``_build_ml_checker``."""

    @patch("deepretro.utils.autosolve.is_valid_smiles", return_value=True)
    def test_keeps_non_hallucinated_pathways(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.return_value = {
            "is_hallucination": False, "probability": 0.1,
        }
        checker = _build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["CC", "O"], ["CCO"]])
        assert status == 200
        assert len(kept) == 2

    @patch("deepretro.utils.autosolve.is_valid_smiles", return_value=True)
    def test_drops_hallucinated_pathways(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.side_effect = [
            {"is_hallucination": True, "probability": 0.9},
            {"is_hallucination": False, "probability": 0.1},
        ]
        checker = _build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["CC"], ["CCO"]])
        assert status == 200
        assert kept == [["CCO"]]

    @patch("deepretro.utils.autosolve.is_valid_smiles", return_value=True)
    def test_joins_reactants_with_dot(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.return_value = {
            "is_hallucination": False, "probability": 0.1,
        }
        checker = _build_ml_checker(mock_clf)
        checker(BENZENE, [["CC", "O"]])
        mock_clf.predict_single.assert_called_once_with(BENZENE, "CC.O")

    @patch("deepretro.utils.autosolve.is_valid_smiles", return_value=False)
    def test_skips_invalid_smiles(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        checker = _build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["INVALID"]])
        assert status == 200
        assert kept == []
        mock_clf.predict_single.assert_not_called()

    @patch("deepretro.utils.autosolve.is_valid_smiles", return_value=True)
    def test_all_hallucinated_returns_empty(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.return_value = {
            "is_hallucination": True, "probability": 0.95,
        }
        checker = _build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["CC"], ["O"]])
        assert status == 200
        assert kept == []


class TestHallucinationModeIntegration:
    """Tests that hallucination_mode correctly wires into the pipeline call."""

    def test_heuristic_passes_true_and_no_fn(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        autosolve(BENZENE, hallucination_mode="heuristic")
        call_kwargs = mock_run.call_args[1]
        assert call_kwargs["hallucination_check"] == "True"
        assert call_kwargs["hallucination_checker_fn"] is None

    def test_none_passes_false_and_no_fn(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        autosolve(BENZENE, hallucination_mode="none")
        call_kwargs = mock_run.call_args[1]
        assert call_kwargs["hallucination_check"] == "False"
        assert call_kwargs["hallucination_checker_fn"] is None

    def test_ml_passes_true_and_a_callable(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        mock_clf = MagicMock()
        with patch(
            "deepretro.utils.autosolve._load_classifier", return_value=mock_clf
        ):
            autosolve(
                BENZENE,
                hallucination_mode="ml",
                hallucination_classifier=mock_clf,
            )
        call_kwargs = mock_run.call_args[1]
        assert call_kwargs["hallucination_check"] == "True"
        assert callable(call_kwargs["hallucination_checker_fn"])

    def test_invalid_mode_raises_before_pipeline(self, mock_pipeline) -> None:
        mock_run = mock_pipeline
        with pytest.raises(ValueError):
            autosolve(BENZENE, hallucination_mode="bad")
        mock_run.assert_not_called()
