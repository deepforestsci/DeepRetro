"""Unit tests for find_optimal_threshold."""

import numpy as np

from deepretro.utils.metrics import find_optimal_threshold


def test_chosen_threshold_maximises_f1():
    """Brute-force verify no other threshold gives a higher F1."""
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
    proba = np.array([0.1] * 85 + [0.4] * 5 + [0.3] * 3 + [0.45] * 4 + [0.9] * 3)
    thr, _ = find_optimal_threshold(y, proba)
    assert thr < 0.5


def test_all_same_probability_does_not_crash():
    """Degenerate input — all probabilities identical — should not error."""
    y = np.array([0, 0, 1, 1])
    proba = np.array([0.5, 0.5, 0.5, 0.5])
    thr, f1 = find_optimal_threshold(y, proba)
    assert isinstance(thr, float)
    assert isinstance(f1, float)
