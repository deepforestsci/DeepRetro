"""Threshold optimization utilities for binary classification."""

from typing import Sequence

import numpy as np


def find_optimal_threshold(
    y_true: Sequence[float],
    probabilities: Sequence[float],
) -> tuple[float, float]:
    """
    Find the classification threshold that maximises F1-score.

    Sweeps the precision-recall curve and picks the threshold where
    the harmonic mean of precision and recall is highest.

    Parameters
    ----------
    y_true : array-like, shape (n_samples,)
        True binary labels (0 or 1).
    probabilities : array-like, shape (n_samples,)
        Predicted probabilities for the positive class.

    Returns
    -------
    threshold : float
        Optimal classification threshold.
    f1 : float
        F1-score at the optimal threshold.

    Examples
    --------
    >>> import numpy as np
    >>> from deepretro.utils.metrics import find_optimal_threshold
    >>> y = np.array([0, 0, 1, 1])
    >>> proba = np.array([0.1, 0.4, 0.6, 0.9])
    >>> thr, f1 = find_optimal_threshold(y, proba)
    >>> 0.0 < thr < 1.0
    True
    >>> f1 > 0.0
    True
    """
    # Import lazily so docs and lightweight tooling can import this module
    # without pulling in sklearn's full scipy stack at module import time.
    from sklearn.metrics import precision_recall_curve

    precision, recall, thresholds = precision_recall_curve(y_true, probabilities)
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
    best_idx = np.argmax(f1_scores)
    if best_idx >= len(thresholds):
        best_idx = len(thresholds) - 1
    return float(thresholds[best_idx]), float(f1_scores[best_idx])
