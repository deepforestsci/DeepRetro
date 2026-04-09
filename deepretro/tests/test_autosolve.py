"""Unit tests for deepretro.algorithms.autosolve.

Tests the autosolve wrapper using mocked pipeline calls so the tests
run without AiZynthFinder models, LLM API keys, or any other heavy
infrastructure.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from deepretro.algorithms.autosolve import autosolve
from deepretro.models.hallucination_helpers import (
    build_ml_checker,
    resolve_hallucination_args,
)

BENZENE = "c1ccccc1"
ETHANOL = "CCO"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

FAKE_TREE = {"type": "mol", "smiles": BENZENE, "children": []}
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
    """Patch rec_run and format_output so tests never touch real infra."""
    with patch(
        "deepretro.algorithms.autosolve.rec_run",
        return_value=(FAKE_TREE, True),
    ) as mock_rec, patch(
        "deepretro.utils.parse.format_output",
        return_value=FAKE_RESULT,
    ) as mock_fmt:
        yield mock_rec, mock_fmt


class TestAutosolve:
    """Tests for ``autosolve``."""

    def test_returns_result_dict(self, mock_pipeline) -> None:
        result = autosolve(BENZENE)
        assert result == FAKE_RESULT

    def test_passes_params_to_rec_run(self, mock_pipeline) -> None:
        mock_rec, _ = mock_pipeline
        autosolve(
            ASPIRIN,
            llm="test-model",
            az_model="USPTO",
            stability_check=False,
            hallucination_mode="none",
            use_protecting_group_feature=True,
        )
        call_kwargs = mock_rec.call_args[1]
        assert call_kwargs["llm"] == "test-model"
        assert call_kwargs["az_model"] == "USPTO"
        assert call_kwargs["stability_flag"] == "False"
        assert call_kwargs["hallucination_check"] == "False"
        assert call_kwargs["use_protecting_group_feature"] is True
        assert call_kwargs["hallucination_checker_fn"] is None

    def test_does_not_write_to_disk(self, mock_pipeline, tmp_path: Path) -> None:
        autosolve(BENZENE)
        assert list(tmp_path.iterdir()) == []


class TestHallucinationModeIntegration:
    """Tests that hallucination_mode correctly wires into the pipeline call."""

    def test_heuristic_passes_true_and_no_fn(self, mock_pipeline) -> None:
        mock_rec, _ = mock_pipeline
        autosolve(BENZENE, hallucination_mode="heuristic")
        call_kwargs = mock_rec.call_args[1]
        assert call_kwargs["hallucination_check"] == "True"
        assert call_kwargs["hallucination_checker_fn"] is None

    def test_none_passes_false_and_no_fn(self, mock_pipeline) -> None:
        mock_rec, _ = mock_pipeline
        autosolve(BENZENE, hallucination_mode="none")
        call_kwargs = mock_rec.call_args[1]
        assert call_kwargs["hallucination_check"] == "False"
        assert call_kwargs["hallucination_checker_fn"] is None

    def test_ml_passes_true_and_a_callable(self, mock_pipeline) -> None:
        mock_rec, _ = mock_pipeline
        mock_clf = MagicMock()
        mock_clf.predict_single = MagicMock()
        with patch(
            "deepretro.algorithms.autosolve.HallucinationClassifier",
            new=type(mock_clf),
        ):
            autosolve(
                BENZENE,
                hallucination_mode="ml",
                hallucination_classifier=mock_clf,
            )
        call_kwargs = mock_rec.call_args[1]
        assert call_kwargs["hallucination_check"] == "True"
        assert callable(call_kwargs["hallucination_checker_fn"])

    def test_invalid_mode_raises(self, mock_pipeline) -> None:
        with pytest.raises(ValueError):
            autosolve(BENZENE, hallucination_mode="bad")


class TestResolveHallucinationArgs:
    """Tests for the standalone resolve_hallucination_args helper."""

    def test_heuristic_mode(self) -> None:
        h_check, h_fn = resolve_hallucination_args("heuristic", None)
        assert h_check == "True"
        assert h_fn is None

    def test_none_mode(self) -> None:
        h_check, h_fn = resolve_hallucination_args("none", None)
        assert h_check == "False"
        assert h_fn is None

    def test_invalid_mode_raises(self) -> None:
        with pytest.raises(ValueError, match="hallucination_mode must be"):
            resolve_hallucination_args("invalid", None)

    def test_ml_mode_without_classifier_raises(self) -> None:
        with pytest.raises(ValueError, match="requires hallucination_classifier"):
            resolve_hallucination_args("ml", None)


class TestBuildMlChecker:
    """Tests for ``build_ml_checker``."""

    @patch("deepretro.utils.utils_molecule.is_valid_smiles", return_value=True)
    def test_keeps_non_hallucinated_pathways(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.return_value = {
            "is_hallucination": False, "probability": 0.1,
        }
        checker = build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["CC", "O"], ["CCO"]])
        assert status == 200
        assert len(kept) == 2

    @patch("deepretro.utils.utils_molecule.is_valid_smiles", return_value=True)
    def test_drops_hallucinated_pathways(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.side_effect = [
            {"is_hallucination": True, "probability": 0.9},
            {"is_hallucination": False, "probability": 0.1},
        ]
        checker = build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["CC"], ["CCO"]])
        assert status == 200
        assert kept == [["CCO"]]

    @patch("deepretro.utils.utils_molecule.is_valid_smiles", return_value=True)
    def test_all_hallucinated_returns_empty(self, _mock_valid) -> None:
        mock_clf = MagicMock()
        mock_clf.predict_single.return_value = {
            "is_hallucination": True, "probability": 0.95,
        }
        checker = build_ml_checker(mock_clf)
        status, kept = checker(BENZENE, [["CC"], ["O"]])
        assert status == 200
        assert kept == []
