"""Unit tests for hallucination helper wiring."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from deepchem.data import NumpyDataset

from deepretro.models.hallucination_helpers import MLChecker, resolve_hallucination


class FixedProbabilityClassifier:
    """Deterministic classifier double for hallucination helper tests."""

    def __init__(self, probability: float) -> None:
        self.threshold = 0.5
        self.probability = probability
        self.featurizer = SimpleNamespace(
            featurize=lambda reactions: np.zeros((len(reactions), 1))
        )

    def predict_probability(
        self,
        dataset: object,
        threshold: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        del threshold
        del dataset
        label = int(self.probability >= self.threshold)
        return np.array([label]), np.array([self.probability])


def test_resolve_hallucination_none_returns_none() -> None:
    """Mode `none` should disable hallucination checking."""
    assert resolve_hallucination("none", None) is None


def test_resolve_hallucination_ml_wraps_classifier_instance() -> None:
    """ML mode should wrap classifier-like instances in MLChecker."""
    checker = resolve_hallucination("ml", FixedProbabilityClassifier(0.1))
    assert isinstance(checker, MLChecker)


def test_resolve_hallucination_heuristic_returns_callable() -> None:
    """Heuristic mode should return the package heuristic checker callable."""
    checker = resolve_hallucination("heuristic", None)
    assert callable(checker)


def test_resolve_hallucination_invalid_mode_raises() -> None:
    """Unsupported modes should raise a clear ValueError."""
    with pytest.raises(ValueError, match="hallucination_mode must be"):
        resolve_hallucination("bad-mode", None)


def test_resolve_hallucination_ml_rejects_invalid_classifier() -> None:
    """ML mode should reject classifier values that are neither path nor model."""
    with pytest.raises(ValueError, match="requires a HallucinationClassifier"):
        resolve_hallucination("ml", object())


def test_mlchecker_keeps_valid_list_pathways() -> None:
    """MLChecker should keep valid list-based pathways with low-risk scores."""
    checker = MLChecker(FixedProbabilityClassifier(0.1))

    status, kept = checker("CCO", [["CC", "O"], ["C", "N"]])

    assert status == 200
    assert kept == [["CC", "O"], ["C", "N"]]


def test_mlchecker_normalizes_valid_string_pathways() -> None:
    """MLChecker should normalize string pathways to one-item pathway lists."""
    checker = MLChecker(FixedProbabilityClassifier(0.1))

    status, kept = checker("CCO", ["CC.O"])

    assert status == 200
    assert kept == [["CC.O"]]


def test_mlchecker_skips_invalid_smiles_and_keeps_remaining_valid_pathways() -> None:
    """MLChecker should skip invalid pathways without failing valid ones."""
    checker = MLChecker(FixedProbabilityClassifier(0.1))

    status, kept = checker("CCO", ["not_a_smiles", ["CC", "O"]])

    assert status == 200
    assert kept == [["CC", "O"]]


def test_mlchecker_returns_400_when_all_pathways_are_predicted_hallucinations() -> None:
    """MLChecker should request retry when all pathways are classified as bad."""
    checker = MLChecker(FixedProbabilityClassifier(0.9))

    status, kept = checker("CCO", [["CC", "O"], ["CC", "N"]])

    assert status == 400
    assert kept == []


def test_mlchecker_returns_400_when_all_pathways_are_invalid() -> None:
    """MLChecker should request retry when every candidate pathway is invalid."""
    checker = MLChecker(FixedProbabilityClassifier(0.1))

    status, kept = checker("CCO", ["not_a_smiles", "still_bad"])

    assert status == 400
    assert kept == []


def test_inline_classifier_matches_predict_probability_contract() -> None:
    """The inline classifier should match the expected tuple return contract."""
    dataset = NumpyDataset(X=np.zeros((1, 1)))
    classifier = FixedProbabilityClassifier(0.1)

    labels, probabilities = classifier.predict_probability(dataset)

    assert labels.shape == (1,)
    assert probabilities.shape == (1,)
