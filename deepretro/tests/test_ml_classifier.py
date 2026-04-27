"""Unit tests for HallucinationClassifier.

Build a small NumpyDataset,
fit, evaluate, predict, and test save/reload roundtrip.
"""

import numpy as np
import pytest
from deepchem.data import NumpyDataset
from deepchem.models import GBDTModel
from unittest.mock import patch

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


# inheritance


def test_classifier_is_gbdtmodel_subclass():
    """HallucinationClassifier must inherit from GBDTModel."""
    assert issubclass(HallucinationClassifier, GBDTModel)


# model_dir paths


def test_user_specified_model_dir(tmp_path):
    """When model_dir is given, model files are saved there."""
    model_dir = str(tmp_path / "user_chosen")
    clf = HallucinationClassifier(model_dir=model_dir, n_estimators=5)
    assert clf.model_dir == model_dir


def test_no_model_dir_still_constructs():
    """When model_dir is None, classifier still constructs without error."""
    clf = HallucinationClassifier(n_estimators=5)
    assert clf.model_dir is not None


def test_load_missing_model_dir(tmp_path):
    """load() from a dir with no model.joblib raises an error."""
    clf = HallucinationClassifier(model_dir=str(tmp_path), n_estimators=5)
    with pytest.raises(Exception):
        clf.load(str(tmp_path))


# evaluate


def test_evaluate_keys_and_scores(trained_clf, toy_dataset):
    scores = trained_clf.evaluate(toy_dataset)
    assert set(scores.keys()) == {
        "roc_auc",
        "accuracy",
        "f1",
        "optimal_threshold",
        "optimal_f1",
    }
    for key in ("roc_auc", "accuracy", "f1", "optimal_f1"):
        assert 0.0 <= scores[key] <= 1.0
    assert 0.0 < scores["optimal_threshold"] < 1.0
    assert trained_clf.threshold == scores["optimal_threshold"]


# predict


def test_predict_probability_shape_and_range(trained_clf, toy_dataset):
    labels, probabilities = trained_clf.predict_probability(toy_dataset)
    assert labels.shape == (len(toy_dataset),)
    assert probabilities.shape == (len(toy_dataset),)
    assert np.all((probabilities >= 0.0) & (probabilities <= 1.0))


def test_predict_probability_override_uses_explicit_threshold():
    clf = HallucinationClassifier(n_estimators=5)
    dataset = NumpyDataset(X=np.zeros((2, 10)))
    with patch.object(
        clf.model, "predict_proba", return_value=np.array([[0.8, 0.2], [0.2, 0.8]])
    ):
        labels, probabilities = clf.predict_probability(
            dataset, threshold=0.7
        )

    np.testing.assert_array_equal(labels, np.array([0, 1]))
    np.testing.assert_array_equal(probabilities, np.array([0.2, 0.8]))


def test_predict_single_override_uses_explicit_threshold():
    clf = HallucinationClassifier(n_estimators=5)
    original_threshold = clf.threshold

    with patch.object(
        clf, "predict_probability", return_value=(np.array([0]), np.array([0.6]))
    ):
        result = clf.predict_single(
            PYRAZOLE_ADDUCT,
            PYRAZOLE_BROMIDE_KETONE,
            threshold=0.7,
        )

    assert result["is_hallucination"] is False
    assert result["probability"] == 0.6
    assert clf.threshold == original_threshold


def test_predict_single_invalid_smiles(trained_clf):
    result = trained_clf.predict_single(INVALID_SMILES, ETHANE_WATER)
    assert result["is_hallucination"] is None
    assert result["probability"] is None
    assert "error" in result


# save / load roundtrip


def test_save_load_roundtrip(trained_clf, toy_dataset, tmp_path):
    trained_clf.evaluate(toy_dataset)
    _, probability_before = trained_clf.predict_probability(toy_dataset)
    saved_threshold = trained_clf.threshold

    save_dir = str(tmp_path / "saved")
    trained_clf.save(save_dir)

    new_clf = HallucinationClassifier(model_dir=save_dir)
    new_clf.load(save_dir)
    assert new_clf.threshold == saved_threshold
    _, probability_after = new_clf.predict_probability(toy_dataset)
    np.testing.assert_array_almost_equal(
        probability_before, probability_after
    )


# reload from saved dir (user flow)


def test_reload_from_saved_dir(trained_clf, toy_dataset, tmp_path):
    """User saves to a dir, then constructs + loads from same dir."""
    save_dir = str(tmp_path / "saved_model")
    trained_clf.evaluate(toy_dataset)
    trained_clf.save(save_dir)

    clf = HallucinationClassifier(model_dir=save_dir)
    clf.load(save_dir)
    assert clf.threshold == trained_clf.threshold

    _, probability_orig = trained_clf.predict_probability(toy_dataset)
    _, probability_loaded = clf.predict_probability(toy_dataset)
    np.testing.assert_array_almost_equal(probability_orig, probability_loaded)
