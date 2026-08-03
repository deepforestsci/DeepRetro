"""Unit tests for :class:`deepretro.models.hallucination_trainer.HallucinationTrainer`.

The constructor, threshold optimisation, config round-trip, and guard-clause
tests are fast (no model fitting). The end-to-end ``train_model`` and
hyperparameter-tuning tests actually fit xgboost models and are marked
``slow``.
"""

import csv
import json
import os

import deepchem as dc
import numpy as np
import pytest
from deepretro.models.hallucination_trainer import HallucinationTrainer

VALID_PRODUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
VALID_REACTANTS = "Cn1nccc1Br.O=C1CCCC[C@H]1O"


# ---------------------------------------------------------------------------
# Constructor branches
# ---------------------------------------------------------------------------


def test_requires_model_type_or_path(tmp_path):
    """Constructing with neither a model type nor a path is an error."""
    with pytest.raises(ValueError, match="model_type.*model_path"):
        HallucinationTrainer(trainer_dir=str(tmp_path))


def test_unknown_model_type_raises(tmp_path):
    """An unregistered model type is rejected."""
    with pytest.raises(ValueError, match="Unknown model_type"):
        HallucinationTrainer(trainer_dir=str(tmp_path), model_type="randomforest")


def test_model_path_without_config_raises(tmp_path):
    """Loading from a directory lacking ``config.json`` fails."""
    with pytest.raises(FileNotFoundError, match="config.json"):
        HallucinationTrainer(trainer_dir=str(tmp_path), model_path=str(tmp_path))


def test_new_xgboost_trainer_populates_config(tmp_path):
    """A fresh xgboost trainer adopts the registry defaults and builds a featurizer."""
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    assert trainer.model_type == "xgboost"
    assert trainer.feat_name == "reactionfeaturizer"
    assert trainer.param_space  # non-empty search space from the registry
    assert trainer.featurizer is not None
    assert os.path.isdir(trainer.model_dir)


def test_model_type_is_lowercased(tmp_path):
    """The model type selector is normalised to lower case."""
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="XGBoost")
    assert trainer.model_type == "xgboost"


# ---------------------------------------------------------------------------
# load_dataset guard
# ---------------------------------------------------------------------------


def test_load_dataset_rejects_multitask(tmp_path):
    """A trainer configured for >1 task refuses the single-label loader."""
    trainer = HallucinationTrainer(
        trainer_dir=str(tmp_path), model_type="xgboost", n_tasks=2
    )
    with pytest.raises(ValueError, match="does not match"):
        trainer.load_dataset("train.csv", "test.csv")


# ---------------------------------------------------------------------------
# evaluate guard
# ---------------------------------------------------------------------------


def test_evaluate_without_model_raises(tmp_path):
    """Evaluating before a model is built is a runtime error."""
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    dummy = dc.data.NumpyDataset(X=np.zeros((2, 3)), y=np.array([[0], [1]]))
    with pytest.raises(RuntimeError, match="not been instantiated"):
        trainer.evaluate(dummy)


# ---------------------------------------------------------------------------
# optimize_threshold (no training required)
# ---------------------------------------------------------------------------


class _FixedProbaModel:
    """Test double whose ``predict`` returns preset positive-class probabilities."""

    def __init__(self, proba):
        # Shape (N, 2): column 1 is the positive-class probability.
        self._pred = np.column_stack([1.0 - np.asarray(proba), np.asarray(proba)])

    def predict(self, dataset):
        return self._pred


def test_optimize_threshold_separates_perfectly(tmp_path):
    """On separable scores the chosen threshold reproduces the true labels."""
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    y = np.array([[0], [0], [1], [1]])
    valid = dc.data.NumpyDataset(X=np.zeros((4, 2)), y=y)
    model = _FixedProbaModel([0.1, 0.2, 0.8, 0.9])

    threshold = trainer.optimize_threshold(model, valid)

    assert isinstance(threshold, float)
    predictions = (np.array([0.1, 0.2, 0.8, 0.9]) >= threshold).astype(int)
    assert list(predictions) == [0, 0, 1, 1]


# ---------------------------------------------------------------------------
# save_config / model_path round-trip (no training required)
# ---------------------------------------------------------------------------


def test_save_config_writes_expected_fields(tmp_path):
    """``save_config`` serialises the trainer state to ``config.json``."""
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    trainer.threshold = 0.42
    trainer.save_config()

    config_path = os.path.join(trainer.model_dir, "config.json")
    with open(config_path) as fh:
        config = json.load(fh)

    assert config["model_type"] == "xgboost"
    assert config["feat_name"] == "reactionfeaturizer"
    assert config["threshold"] == pytest.approx(0.42)
    assert config["n_tasks"] == 1


def test_model_path_constructor_restores_config(tmp_path):
    """A trainer built from a saved ``model_path`` restores its configuration."""
    src = HallucinationTrainer(trainer_dir=str(tmp_path / "src"), model_type="xgboost")
    src.threshold = 0.37
    src.save_config()

    restored = HallucinationTrainer(
        trainer_dir=str(tmp_path / "dst"), model_path=src.model_dir
    )
    assert restored.model_type == "xgboost"
    assert restored.threshold == pytest.approx(0.37)
    assert restored.feat_name == "reactionfeaturizer"


# ---------------------------------------------------------------------------
# parameter_tuning early-exit (no training required)
# ---------------------------------------------------------------------------


def test_parameter_tuning_without_space_returns_defaults(tmp_path):
    """With no search space, tuning short-circuits to the current settings."""
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    trainer.param_space = {}
    dummy = dc.data.NumpyDataset(X=np.zeros((2, 3)), y=np.array([[0], [1]]))

    params, threshold = trainer.parameter_tuning(dummy)

    assert params == trainer.model_params
    assert threshold == trainer.threshold


# ---------------------------------------------------------------------------
# End-to-end training (slow)
# ---------------------------------------------------------------------------


def _write_reaction_csv(path, rows):
    """Write a product/reactants/label CSV the trainer can load."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["product", "reactants", "label"])
        writer.writerows(rows)


@pytest.fixture(scope="module")
def reaction_csvs(tmp_path_factory):
    """Create tiny balanced train/test CSVs for end-to-end training."""
    rows = [
        (VALID_PRODUCT, VALID_REACTANTS, 0),
        ("CCO", "CC=O", 0),
        ("CC(=O)O", "CC(=O)Cl.O", 0),
        ("c1ccccc1O", "c1ccccc1Br.O", 0),
        ("CCN", "CC#N", 1),
        ("CCCBr", "CCCO", 1),
        ("c1ccccc1N", "c1ccccc1", 1),
        ("CC(C)O", "CC(C)Cl", 1),
    ]
    workdir = tmp_path_factory.mktemp("trainer_csvs")
    train_csv = workdir / "train.csv"
    test_csv = workdir / "test.csv"
    _write_reaction_csv(train_csv, rows)
    _write_reaction_csv(test_csv, rows)
    return str(train_csv), str(test_csv)


@pytest.mark.slow
def test_train_model_end_to_end(tmp_path, reaction_csvs):
    """A full train run returns a fitted model, threshold metrics, and a saved config."""
    train_csv, test_csv = reaction_csvs
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    train_ds, test_ds = trainer.load_dataset(train_csv, test_csv)

    model, scores = trainer.train_model(train_ds, test_ds)

    assert model is not None
    for key in ("threshold_accuracy", "threshold_precision", "threshold_recall"):
        assert key in scores
    assert os.path.exists(os.path.join(trainer.model_dir, "config.json"))


@pytest.mark.slow
def test_train_model_with_tuning(tmp_path, reaction_csvs):
    """Training with ``tune_params`` runs the single-split search and updates state."""
    train_csv, test_csv = reaction_csvs
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    train_ds, test_ds = trainer.load_dataset(train_csv, test_csv)

    model, _scores = trainer.train_model(
        train_ds, test_ds, tune_params=True, n_trials=1, k_folds=1
    )

    assert model is not None
    assert 0.0 <= trainer.threshold <= 1.0
