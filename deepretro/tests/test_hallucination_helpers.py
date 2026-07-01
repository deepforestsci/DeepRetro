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


def _drop_second_checker(product: str, pathways: list) -> tuple[int, list]:
    """Order-preserving checker double: keep every pathway except index 1."""
    del product
    kept = [pathway for i, pathway in enumerate(pathways) if i != 1]
    if kept:
        return 200, kept
    return 400, []


class TestFilterWithChecker:
    """Tests for ``hallucination_helpers.filter_with_checker`` realignment."""

    def test_none_checker_passes_everything_through(self) -> None:
        from deepretro.models.hallucination_helpers import filter_with_checker

        pathways = [["CC"], ["O"]]
        result = filter_with_checker(
            "CCO", pathways, ["a", "b"], [0.1, 0.2], None
        )
        assert result == (pathways, ["a", "b"], [0.1, 0.2])

    def test_kept_subset_keeps_aligned_explanations_and_confidence(self) -> None:
        from deepretro.models.hallucination_helpers import filter_with_checker

        pathways = [["A"], ["B"], ["C"]]
        kept, expl, conf = filter_with_checker(
            "P",
            pathways,
            ["expl_a", "expl_b", "expl_c"],
            [0.1, 0.2, 0.3],
            _drop_second_checker,
        )
        assert kept == [["A"], ["C"]]
        assert expl == ["expl_a", "expl_c"]
        assert conf == [0.1, 0.3]

    def test_rejected_all_returns_empty_lists(self) -> None:
        from deepretro.models.hallucination_helpers import filter_with_checker

        def reject_all(product: str, pathways: list) -> tuple[int, list]:
            del product, pathways
            return 400, []

        assert filter_with_checker(
            "P", [["A"]], ["e"], [0.5], reject_all
        ) == ([], [], [])

    def test_length_mismatch_raises(self) -> None:
        from deepretro.models.hallucination_helpers import filter_with_checker

        with pytest.raises(ValueError, match="same length"):
            filter_with_checker("P", [["A"], ["B"]], ["only_one"], [0.1, 0.2], None)

    def test_unknown_returned_pathway_raises(self) -> None:
        from deepretro.models.hallucination_helpers import filter_with_checker

        def invent_pathway(product: str, pathways: list) -> tuple[int, list]:
            del product, pathways
            return 200, [["NOT", "IN", "INPUT"]]

        with pytest.raises(ValueError, match="not present in original"):
            filter_with_checker("P", [["A"]], ["e"], [0.5], invent_pathway)
