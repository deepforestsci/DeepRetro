"""Unit tests for find_optimal_threshold."""

import numpy as np
import pytest

from deepretro.utils.metrics import find_optimal_threshold


def test_return_types_and_bounds():
    y = np.array([0, 0, 1, 1])
    proba = np.array([0.1, 0.4, 0.6, 0.9])
    thr, f1 = find_optimal_threshold(y, proba)
    assert isinstance(thr, float)
    assert isinstance(f1, float)
    assert 0.0 < thr < 1.0
    assert 0.0 < f1 <= 1.0


def test_perfect_predictions_f1_is_1():
    y = np.array([0, 0, 0, 1, 1, 1])
    proba = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    _, f1 = find_optimal_threshold(y, proba)
    assert f1 == pytest.approx(1.0)


def test_random_predictions_low_f1():
    rng = np.random.RandomState(42)
    y = np.array([0] * 50 + [1] * 50)
    proba = rng.rand(100)  # random probabilities, no signal
    _, f1 = find_optimal_threshold(y, proba)
    assert f1 < 0.8


def test_chosen_threshold_maximises_f1():
    """Verify that no other threshold gives a higher F1."""
    rng = np.random.RandomState(99)
    y = np.array([0] * 80 + [1] * 20)
    proba = np.where(y == 1, rng.uniform(0.5, 1.0, 100), rng.uniform(0.0, 0.5, 100))
    thr, best_f1 = find_optimal_threshold(y, proba)

    for t in np.linspace(0.01, 0.99, 200):
        preds = (proba >= t).astype(int)
        tp = ((preds == 1) & (y == 1)).sum()
        fp = ((preds == 1) & (y == 0)).sum()
        fn = ((preds == 0) & (y == 1)).sum()
        if tp + fp == 0 or tp + fn == 0:
            continue
        p = tp / (tp + fp)
        r = tp / (tp + fn)
        f1 = 2 * p * r / (p + r)
        assert f1 <= best_f1 + 1e-6


def test_imbalanced_data_threshold_below_05():
    """With rare positives, optimal threshold should shift below 0.5."""
    y = np.array([0] * 90 + [1] * 10)
    proba = np.array(
        [0.1] * 85 + [0.4] * 5
        + [0.3] * 3 + [0.45] * 4 + [0.9] * 3
    )
    thr, _ = find_optimal_threshold(y, proba)
    assert thr < 0.5


def test_all_same_probability_does_not_crash():
    y = np.array([0, 0, 1, 1])
    proba = np.array([0.5, 0.5, 0.5, 0.5])
    thr, f1 = find_optimal_threshold(y, proba)
    assert isinstance(thr, float)
    assert isinstance(f1, float)
