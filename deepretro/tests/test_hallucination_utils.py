"""Tests for the shared hallucination-model utilities."""

import deepchem as dc
import pytest

from deepretro.models.hallucination_utils import (
    MODEL_CONFIGS,
    KFoldRandomHyperparamOpt,
    create_model_instance,
)


def test_model_configs_registers_xgboost():
    assert "xgboost" in MODEL_CONFIGS


def test_xgboost_config_has_required_keys():
    config = MODEL_CONFIGS["xgboost"]
    for key in (
        "model_class",
        "base_estimator",
        "is_sklearn",
        "default_params",
        "param_space",
        "default_feat",
        "default_feat_params",
    ):
        assert key in config, f"missing key: {key}"
    assert config["is_sklearn"] is True
    assert config["default_feat"] == "reactionfeaturizer"


def test_create_model_instance_returns_sklearn_model(tmp_path):
    model = create_model_instance(
        "xgboost", model_dir=str(tmp_path), n_tasks=1, max_depth=3
    )
    assert isinstance(model, dc.models.SklearnModel)


def test_create_model_instance_sets_model_dir(tmp_path):
    model_dir = str(tmp_path / "xgb")
    model = create_model_instance("xgboost", model_dir=model_dir)
    assert model.model_dir == model_dir


def test_create_model_instance_forwards_params_to_estimator(tmp_path):
    model = create_model_instance("xgboost", model_dir=str(tmp_path), max_depth=7)
    assert model.model.get_params()["max_depth"] == 7


def test_create_model_instance_pops_infra_kwargs(tmp_path):
    model = create_model_instance("xgboost", model_dir=str(tmp_path), n_tasks=1)
    assert "n_tasks" not in model.model.get_params()


def test_create_model_instance_is_case_insensitive(tmp_path):
    model = create_model_instance("XGBoost", model_dir=str(tmp_path))
    assert isinstance(model, dc.models.SklearnModel)


def test_create_model_instance_unknown_type_raises(tmp_path):
    with pytest.raises(KeyError):
        create_model_instance("randomforest", model_dir=str(tmp_path))


def test_kfold_opt_is_random_opt_subclass():
    assert issubclass(KFoldRandomHyperparamOpt, dc.hyper.RandomHyperparamOpt)


def test_hyperparam_kfold_search_preserves_selected_params_and_score_key(
    tmp_path, monkeypatch
):
    class FixedScoreModel:
        def __init__(self, params, model_score):
            self.params = params
            self.model_score = model_score

        def fit(self, dataset, **kwargs):
            pass

        def evaluate(self, dataset, metrics):
            return {metrics[0].name: self.model_score}

    fold_scores = iter([0.6, 0.9])

    def builder(**params):
        return FixedScoreModel(params.copy(), next(fold_scores))

    optimizer = KFoldRandomHyperparamOpt(model_builder=builder, max_iter=1)
    selected_params = {"max_depth": 3, "n_tasks": 1}
    monkeypatch.setattr(
        optimizer,
        "generate_random_hyperparam_values",
        lambda params_dict, max_iter: [selected_params.copy()],
    )
    metric = dc.metrics.Metric(dc.metrics.roc_auc_score)

    best_models, best_hyperparams, all_scores = optimizer.hyperparam_kfold_search(
        params_dict={"max_depth": [3], "n_tasks": [1]},
        folds=[(object(), object()), (object(), object())],
        metric=metric,
        logdir=str(tmp_path),
    )

    score_key = dc.hyper.base_classes._convert_hyperparam_dict_to_filename(
        selected_params
    )
    assert best_hyperparams == selected_params
    assert "model_dir" not in best_hyperparams
    assert all_scores == {score_key: pytest.approx(0.75)}
    assert all_scores[
        dc.hyper.base_classes._convert_hyperparam_dict_to_filename(best_hyperparams)
    ] == pytest.approx(0.75)
    assert len(best_models) == 2
    assert best_models[0].params["model_dir"] != best_models[1].params["model_dir"]
