"""Unit tests for HallucinationClassifier.

Build a small NumpyDataset,
fit, evaluate, predict, and test save/reload roundtrip.
"""

import json

import numpy as np
import pytest
from deepchem.data import NumpyDataset

from deepretro.models.hallucination_classifier import HallucinationClassifier

PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"
ETHANOL = "CCO"
ETHANE_WATER = "CC.O"
INVALID_SMILES = "not_a_smiles!!!"


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


# __init__

def test_default_threshold():
    clf = HallucinationClassifier()
    assert clf.threshold == 0.5


# fit + evaluate

def test_evaluate_returns_expected_keys(trained_clf, toy_dataset):
    scores = trained_clf.evaluate(toy_dataset)
    assert set(scores.keys()) == {
        "roc_auc", "accuracy", "f1", "optimal_threshold", "optimal_f1",
    }


def test_evaluate_scores_are_valid(trained_clf, toy_dataset):
    scores = trained_clf.evaluate(toy_dataset)
    for key in ("roc_auc", "accuracy", "f1", "optimal_f1"):
        assert 0.0 <= scores[key] <= 1.0
    assert 0.0 < scores["optimal_threshold"] < 1.0


def test_evaluate_updates_threshold(trained_clf, toy_dataset):
    scores = trained_clf.evaluate(toy_dataset)
    assert trained_clf.threshold == scores["optimal_threshold"]


# predict_proba + predict

def test_predict_proba_shape_and_range(trained_clf, toy_dataset):
    proba = trained_clf.predict_proba(toy_dataset)
    assert proba.shape == (len(toy_dataset),)
    assert np.all((proba >= 0.0) & (proba <= 1.0))


def test_predict_returns_binary_labels(trained_clf, toy_dataset):
    labels, proba = trained_clf.predict(toy_dataset)
    assert set(np.unique(labels)).issubset({0, 1})
    assert labels.shape == proba.shape


# predict_single

def test_predict_single_invalid_smiles(trained_clf):
    result = trained_clf.predict_single(INVALID_SMILES, ETHANE_WATER)
    assert result["is_hallucination"] is None
    assert result["probability"] is None
    assert "error" in result


# save / load roundtrip

def test_save_load_preserves_threshold(trained_clf, toy_dataset, tmp_path):
    # evaluate to set a non-default threshold
    trained_clf.evaluate(toy_dataset)
    saved_threshold = trained_clf.threshold

    # save
    save_dir = str(tmp_path / "saved")
    trained_clf.save(save_dir)

    # metadata.json written
    with open(tmp_path / "saved" / "metadata.json") as f:
        meta = json.load(f)
    assert meta["optimal_threshold"] == saved_threshold

    # reload into a fresh classifier
    new_clf = HallucinationClassifier(model_dir=save_dir)
    new_clf.load(save_dir)
    assert new_clf.threshold == saved_threshold


def test_save_load_predictions_match(trained_clf, toy_dataset, tmp_path):
    proba_before = trained_clf.predict_proba(toy_dataset)

    save_dir = str(tmp_path / "saved2")
    trained_clf.save(save_dir)

    new_clf = HallucinationClassifier(model_dir=save_dir)
    new_clf.load(save_dir)
    proba_after = new_clf.predict_proba(toy_dataset)

    np.testing.assert_array_almost_equal(proba_before, proba_after)


# from_pretrained (uses local cache, skips S3 download)

def test_from_pretrained_loads_from_cache(trained_clf, toy_dataset, tmp_path):
    """Save model to a dir, then from_pretrained should load from cache (no download)."""
    cache_dir = tmp_path / "pretrained_cache"
    trained_clf.save(str(cache_dir))

    # files already in cache → urlretrieve never called
    clf = HallucinationClassifier.from_pretrained(cache_dir=str(cache_dir))
    assert clf.threshold == trained_clf.threshold

    proba_orig = trained_clf.predict_proba(toy_dataset)
    proba_loaded = clf.predict_proba(toy_dataset)
    np.testing.assert_array_almost_equal(proba_orig, proba_loaded)

