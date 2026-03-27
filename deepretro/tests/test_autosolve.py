"""Unit tests for deepretro.utils.autosolve.

Tests the autosolve and autosolve_async wrappers using mocked pipeline
calls so the tests run without AiZynthFinder models, LLM API keys, or
any other heavy infrastructure.
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from deepretro.utils.autosolve import _safe_filename, autosolve, autosolve_async

BENZENE = "c1ccccc1"
ETHANOL = "CCO"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

FAKE_RESULT = {
    "dependencies": {"1": []},
    "steps": [
        {
            "step": "1",
            "reactants": [{"smiles": ETHANOL}],
            "products": [{"smiles": BENZENE}],
        }
    ],
}


class TestSafeFilename:
    """Tests for ``_safe_filename``."""

    def test_simple_smiles(self) -> None:
        name = _safe_filename(ETHANOL)
        assert ETHANOL in name
        assert len(name) <= 100

    def test_slashes_replaced(self) -> None:
        smiles = "C/C=C/C"
        name = _safe_filename(smiles)
        assert "/" not in name
        assert "\\" not in name

    def test_long_smiles_truncated(self) -> None:
        long_smiles = "C" * 200
        name = _safe_filename(long_smiles)
        assert len(name) <= 100

    def test_same_input_same_output(self) -> None:
        assert _safe_filename(BENZENE) == _safe_filename(BENZENE)

    def test_different_input_different_output(self) -> None:
        assert _safe_filename(BENZENE) != _safe_filename(ETHANOL)


class TestAutosolve:
    """Tests for ``autosolve``."""

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_returns_result_dict(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        result = autosolve(BENZENE)
        assert result == FAKE_RESULT
        mock_run.assert_called_once()

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_clears_cache_by_default(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        autosolve(BENZENE)
        mock_clear.assert_called_once_with(BENZENE)

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_skips_cache_clear_when_disabled(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        autosolve(BENZENE, clear_cache=False)
        mock_clear.assert_not_called()

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_passes_all_params_to_pipeline(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        autosolve(
            ASPIRIN,
            llm="test-model",
            az_model="USPTO",
            stability_flag="False",
            hallucination_check="False",
            use_protecting_group_feature=True,
        )
        mock_run.assert_called_once_with(
            ASPIRIN,
            llm="test-model",
            az_model="USPTO",
            stability_flag="False",
            hallucination_check="False",
            use_protecting_group_feature=True,
        )

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_writes_json_when_output_path_given(
        self, mock_clear: MagicMock, mock_run: MagicMock, tmp_path: Path
    ) -> None:
        out = tmp_path / "result.json"
        autosolve(BENZENE, output_path=out)
        assert out.exists()
        data = json.loads(out.read_text())
        assert data == FAKE_RESULT

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_creates_parent_dirs_for_output(
        self, mock_clear: MagicMock, mock_run: MagicMock, tmp_path: Path
    ) -> None:
        out = tmp_path / "nested" / "deep" / "result.json"
        autosolve(BENZENE, output_path=out)
        assert out.exists()

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_no_file_written_when_output_path_is_none(
        self, mock_clear: MagicMock, mock_run: MagicMock, tmp_path: Path
    ) -> None:
        autosolve(BENZENE)
        assert list(tmp_path.iterdir()) == []


class TestAutosolveAsync:
    """Tests for ``autosolve_async``."""

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_returns_one_result_per_input(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        results = asyncio.run(autosolve_async([BENZENE, ETHANOL]))
        assert len(results) == 2
        assert all(r == FAKE_RESULT for r in results)

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_empty_list_returns_empty(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        results = asyncio.run(autosolve_async([]))
        assert results == []

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_clears_cache_for_each_molecule(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        asyncio.run(autosolve_async([BENZENE, ETHANOL, ASPIRIN]))
        assert mock_clear.call_count == 3

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_skips_cache_clear_when_disabled(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        asyncio.run(autosolve_async([BENZENE], clear_cache=False))
        mock_clear.assert_not_called()

    @patch("deepretro.utils.autosolve._run_retrosynthesis", side_effect=RuntimeError("LLM down"))
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_failed_molecule_returns_error_dict(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        results = asyncio.run(autosolve_async([BENZENE]))
        assert len(results) == 1
        assert "error" in results[0]
        assert results[0]["smiles"] == BENZENE
        assert "LLM down" in results[0]["error"]

    @patch("deepretro.utils.autosolve._run_retrosynthesis")
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_partial_failure_does_not_crash_batch(
        self, mock_clear: MagicMock, mock_run: MagicMock
    ) -> None:
        mock_run.side_effect = [FAKE_RESULT, RuntimeError("fail"), FAKE_RESULT]
        results = asyncio.run(autosolve_async([BENZENE, ETHANOL, ASPIRIN]))
        assert len(results) == 3
        assert results[0] == FAKE_RESULT
        assert "error" in results[1]
        assert results[2] == FAKE_RESULT

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_writes_output_files_when_dir_given(
        self, mock_clear: MagicMock, mock_run: MagicMock, tmp_path: Path
    ) -> None:
        asyncio.run(autosolve_async([BENZENE, ETHANOL], output_dir=tmp_path))
        json_files = list(tmp_path.glob("*.json"))
        assert len(json_files) == 2
        for f in json_files:
            data = json.loads(f.read_text())
            assert data == FAKE_RESULT

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_no_files_written_when_dir_is_none(
        self, mock_clear: MagicMock, mock_run: MagicMock, tmp_path: Path
    ) -> None:
        asyncio.run(autosolve_async([BENZENE]))
        assert list(tmp_path.iterdir()) == []

    @patch("deepretro.utils.autosolve._run_retrosynthesis", return_value=FAKE_RESULT)
    @patch("deepretro.utils.autosolve.clear_cache_for_molecule")
    def test_respects_max_concurrent(self, mock_clear: MagicMock, mock_run: MagicMock) -> None:
        molecules = [f"C{'C' * i}O" for i in range(10)]
        results = asyncio.run(autosolve_async(molecules, max_concurrent=2))
        assert len(results) == 10
        assert mock_run.call_count == 10
