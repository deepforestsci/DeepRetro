"""Unit tests for HallucinationClassifier.

Build a small NumpyDataset,
fit, evaluate, predict, and test save/reload roundtrip.
"""

import numpy as np
import pytest
from deepchem.data import NumpyDataset

from deepretro.models.hallucination_classifier import HallucinationClassifier

PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"
INVALID_SMILES = "not_a_smiles!!!"
ETHANE_WATER = "CC.O"


@pytest.fixture(scope="module")
def toy_dataset():
    """Small synthetic dataset, enough for GBDTModel's 80/20 internal split."""
    rng = np.random.RandomState(42)
    X = rng.rand(40, 10)
    y = np.array([0] * 20 + [1] * 20).reshape(-1, 1)
    return NumpyDataset(X=X, y=y)


@pytest.fixture(scope="module")
def trained_clf(toy_dataset, tmp_path_factory):
    """Train once, reuse across all tests in this module."""
    model_dir = str(tmp_path_factory.mktemp("model"))
    clf = HallucinationClassifier(
        model_dir=model_dir,
        n_estimators=10,
        early_stopping_rounds=5,
    )
    clf.fit(toy_dataset)
    return clf


# evaluate


def test_evaluate_keys_and_scores(trained_clf, toy_dataset):
    scores = trained_clf.evaluate(toy_dataset)
    assert set(scores.keys()) == {
        "roc_auc", "accuracy", "f1", "optimal_threshold", "optimal_f1",
    }
    for key in ("roc_auc", "accuracy", "f1", "optimal_f1"):
        assert 0.0 <= scores[key] <= 1.0
    assert 0.0 < scores["optimal_threshold"] < 1.0
    assert trained_clf.threshold == scores["optimal_threshold"]


# predict


def test_predict_probability_shape_and_range(trained_clf, toy_dataset):
    probability = trained_clf.predict_probability(toy_dataset)
    assert probability.shape == (len(toy_dataset),)
    assert np.all((probability >= 0.0) & (probability <= 1.0))


def test_predict_single_invalid_smiles(trained_clf):
    result = trained_clf.predict_single(INVALID_SMILES, ETHANE_WATER)
    assert result["is_hallucination"] is None
    assert result["probability"] is None
    assert "error" in result


# save / load roundtrip


def test_save_load_roundtrip(trained_clf, toy_dataset, tmp_path):
    trained_clf.evaluate(toy_dataset)
    probability_before = trained_clf.predict_probability(toy_dataset)
    saved_threshold = trained_clf.threshold

    save_dir = str(tmp_path / "saved")
    trained_clf.save(save_dir)

    new_clf = HallucinationClassifier(model_dir=save_dir)
    new_clf.load(save_dir)
    assert new_clf.threshold == saved_threshold
    np.testing.assert_array_almost_equal(
        probability_before, new_clf.predict_probability(toy_dataset)
    )


# from_pretrained


def test_from_pretrained_loads_from_local_dir(trained_clf, toy_dataset, tmp_path):
    """Save model to a dir, then from_pretrained loads it back."""
    save_dir = tmp_path / "pretrained_local"
    trained_clf.save(str(save_dir))

    clf = HallucinationClassifier.from_pretrained(model_dir=str(save_dir))
    assert clf.threshold == trained_clf.threshold

    probability_orig = trained_clf.predict_probability(toy_dataset)
    probability_loaded = clf.predict_probability(toy_dataset)
    np.testing.assert_array_almost_equal(probability_orig, probability_loaded)
