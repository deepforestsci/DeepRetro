"""Tests for :class:`deepretro.models.hallucination_trainer.HallucinationTrainer`."""

import json
import os

import deepchem as dc
import numpy as np
import pytest

import deepretro.models.hallucination_trainer as trainer_module
from deepretro.models.hallucination_trainer import HallucinationTrainer


def test_requires_model_type_or_path(tmp_path):
    with pytest.raises(ValueError, match="model_type.*model_path"):
        HallucinationTrainer(trainer_dir=str(tmp_path))


def test_unknown_model_type_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown model_type"):
        HallucinationTrainer(trainer_dir=str(tmp_path), model_type="randomforest")


def test_model_path_without_config_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="config.json"):
        HallucinationTrainer(trainer_dir=str(tmp_path), model_path=str(tmp_path))


def test_new_xgboost_trainer_populates_config(tmp_path):
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    assert trainer.model_type == "xgboost"
    assert trainer.feat_name == "reactionfeaturizer"
    assert trainer.param_space
    assert trainer.featurizer is not None
    assert os.path.isdir(trainer.model_dir)


def test_model_type_is_lowercased(tmp_path):
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="XGBoost")
    assert trainer.model_type == "xgboost"


def test_load_dataset_rejects_multitask(tmp_path):
    trainer = HallucinationTrainer(
        trainer_dir=str(tmp_path), model_type="xgboost", n_tasks=2
    )
    with pytest.raises(ValueError, match="does not match"):
        trainer.load_dataset("train.csv", "test.csv")


def test_evaluate_without_model_raises(tmp_path):
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    dummy = dc.data.NumpyDataset(X=np.zeros((2, 3)), y=np.array([[0], [1]]))
    with pytest.raises(RuntimeError, match="not been instantiated"):
        trainer.evaluate(dummy)


class _FixedProbaModel:
    def __init__(self, probabilities):
        probabilities = np.asarray(probabilities)
        self._prediction = np.column_stack([1.0 - probabilities, probabilities])

    def predict(self, dataset):
        return self._prediction


def test_optimize_threshold_separates_perfectly(tmp_path):
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    y = np.array([[0], [0], [1], [1]])
    valid = dc.data.NumpyDataset(X=np.zeros((4, 2)), y=y)
    model = _FixedProbaModel([0.1, 0.2, 0.8, 0.9])

    threshold = trainer.optimize_threshold(model, valid)

    assert isinstance(threshold, float)
    predictions = (np.array([0.1, 0.2, 0.8, 0.9]) >= threshold).astype(int)
    assert list(predictions) == [0, 0, 1, 1]


def test_save_config_writes_expected_fields(tmp_path):
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
    src = HallucinationTrainer(trainer_dir=str(tmp_path / "src"), model_type="xgboost")
    src.threshold = 0.37
    src.save_config()

    restored = HallucinationTrainer(
        trainer_dir=str(tmp_path / "dst"), model_path=src.model_dir
    )
    assert restored.model_type == "xgboost"
    assert restored.threshold == pytest.approx(0.37)
    assert restored.feat_name == "reactionfeaturizer"


def test_parameter_tuning_without_space_returns_defaults(tmp_path):
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    trainer.param_space = {}
    dummy = dc.data.NumpyDataset(X=np.zeros((2, 3)), y=np.array([[0], [1]]))

    params, threshold = trainer.parameter_tuning(dummy)

    assert params == trainer.model_params
    assert threshold == trainer.threshold


def test_parameter_tuning_reports_score_for_selected_kfold_params(
    tmp_path, monkeypatch, capsys
):
    trainer = HallucinationTrainer(trainer_dir=str(tmp_path), model_type="xgboost")
    selected_params = {"max_depth": 7, "n_tasks": 1}
    score_key = dc.hyper.base_classes._convert_hyperparam_dict_to_filename(
        selected_params
    )

    class FixedSplitter:
        def k_fold_split(self, dataset, k_folds):
            return [(object(), object()) for _ in range(k_folds)]

    class FixedOptimizer:
        def __init__(self, model_builder, max_iter):
            pass

        def hyperparam_kfold_search(self, **kwargs):
            polluted_params = selected_params | {"model_dir": "final-fold"}
            return [object(), object()], polluted_params, {score_key: 0.875}

    monkeypatch.setattr(
        trainer_module.dc.splits,
        "RandomStratifiedSplitter",
        lambda: FixedSplitter(),
    )
    monkeypatch.setattr(trainer_module, "KFoldRandomHyperparamOpt", FixedOptimizer)
    monkeypatch.setattr(trainer, "optimize_threshold", lambda model, dataset: 0.4)

    best_params, threshold = trainer.parameter_tuning(object(), n_trials=1, k_folds=2)

    assert best_params == {"max_depth": 7}
    assert threshold == pytest.approx(0.4)
    assert "PRC-AUC Score: 0.8750" in capsys.readouterr().out
