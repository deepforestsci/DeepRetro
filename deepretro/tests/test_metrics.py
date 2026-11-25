"""Unit tests for find_optimal_threshold."""

import numpy as np
import pytest

from deepretro.utils.metrics import find_optimal_threshold


def test_chosen_threshold_maximises_f1() -> None:
    """No threshold on a fine grid gives a higher F1."""
    rng = np.random.RandomState(99)
    y = np.array([0] * 80 + [1] * 20)
    proba = np.where(y == 1, rng.uniform(0.5, 1.0, 100), rng.uniform(0.0, 0.5, 100))
    thr, best_f1 = find_optimal_threshold(y.tolist(), proba.tolist())

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


def test_rare_positives_pick_threshold_below_half() -> None:
    y = [0] * 90 + [1] * 10
    proba = [0.1] * 85 + [0.4] * 5 + [0.3] * 3 + [0.45] * 4 + [0.9] * 3
    thr, f1 = find_optimal_threshold(y, proba)
    assert thr == 0.45
    # 7 of 10 positives score >= 0.45, no false positives
    assert f1 == pytest.approx(14 / 17)


def test_identical_probabilities_call_every_sample_positive() -> None:
    thr, f1 = find_optimal_threshold([0, 0, 1, 1], [0.5, 0.5, 0.5, 0.5])
    assert thr == 0.5
    assert f1 == pytest.approx(2 / 3)  # precision 1/2, recall 1


def test_tie_selects_lowest_threshold() -> None:
    """When several thresholds give the same best F1, the lowest one wins."""
    thr, f1 = find_optimal_threshold([1, 0, 0, 1], [0.2, 0.4, 0.6, 0.8])
    assert thr == 0.2
    assert f1 == pytest.approx(2 / 3)
