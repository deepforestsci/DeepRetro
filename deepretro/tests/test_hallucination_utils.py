"""Unit tests for :mod:`deepretro.models.hallucination_utils`.

Covers the model-configuration registry, the ``create_model_instance``
factory, and the K-Fold hyperparameter optimiser. These are the shared
building blocks used by :class:`~deepretro.models.hallucination_trainer.HallucinationTrainer`
and :class:`~deepretro.models.hallucination_checker.HallucinationChecker`.
"""

import deepchem as dc
import numpy as np
import pytest
from deepretro.models.hallucination_utils import (
    MODEL_CONFIGS,
    KFoldRandomHyperparamOpt,
    create_model_instance,
)

# ---------------------------------------------------------------------------
# MODEL_CONFIGS registry
# ---------------------------------------------------------------------------


def test_model_configs_registers_xgboost():
    """The registry must expose an ``xgboost`` configuration."""
    assert "xgboost" in MODEL_CONFIGS


def test_xgboost_config_has_required_keys():
    """The xgboost config must carry every key the factory/trainer read."""
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


# ---------------------------------------------------------------------------
# create_model_instance
# ---------------------------------------------------------------------------


def test_create_model_instance_returns_sklearn_model(tmp_path):
    """An xgboost request returns a DeepChem ``SklearnModel`` wrapper."""
    model = create_model_instance(
        "xgboost", model_dir=str(tmp_path), n_tasks=1, max_depth=3
    )
    assert isinstance(model, dc.models.SklearnModel)


def test_create_model_instance_sets_model_dir(tmp_path):
    """The provided ``model_dir`` is threaded through to the wrapper."""
    model_dir = str(tmp_path / "xgb")
    model = create_model_instance("xgboost", model_dir=model_dir)
    assert model.model_dir == model_dir


def test_create_model_instance_forwards_params_to_estimator(tmp_path):
    """Extra kwargs are forwarded to the wrapped XGBClassifier."""
    model = create_model_instance("xgboost", model_dir=str(tmp_path), max_depth=7)
    assert model.model.get_params()["max_depth"] == 7


def test_create_model_instance_pops_infra_kwargs(tmp_path):
    """``n_tasks`` is infrastructure metadata and must not reach the estimator."""
    model = create_model_instance("xgboost", model_dir=str(tmp_path), n_tasks=1)
    assert "n_tasks" not in model.model.get_params()


def test_create_model_instance_is_case_insensitive(tmp_path):
    """The model type lookup lower-cases its argument."""
    model = create_model_instance("XGBoost", model_dir=str(tmp_path))
    assert isinstance(model, dc.models.SklearnModel)


def test_create_model_instance_unknown_type_raises(tmp_path):
    """An unregistered model type raises ``KeyError``."""
    with pytest.raises(KeyError):
        create_model_instance("randomforest", model_dir=str(tmp_path))


# ---------------------------------------------------------------------------
# KFoldRandomHyperparamOpt
# ---------------------------------------------------------------------------


def _balanced_dataset(
    n_per_class: int, n_features: int, seed: int
) -> dc.data.NumpyDataset:
    """Build a tiny balanced binary NumpyDataset for fast xgboost fitting."""
    rng = np.random.RandomState(seed)
    x_neg = rng.rand(n_per_class, n_features)
    x_pos = rng.rand(n_per_class, n_features) + 1.0
    x = np.vstack([x_neg, x_pos])
    y = np.array([0] * n_per_class + [1] * n_per_class).reshape(-1, 1)
    return dc.data.NumpyDataset(X=x, y=y)


def test_kfold_opt_is_random_opt_subclass():
    """The optimiser extends DeepChem's ``RandomHyperparamOpt``."""
    assert issubclass(KFoldRandomHyperparamOpt, dc.hyper.RandomHyperparamOpt)


@pytest.mark.slow
def test_hyperparam_kfold_search_returns_models_per_fold(tmp_path):
    """The K-Fold search returns one trained model per fold plus scores."""
    folds = [
        (
            _balanced_dataset(6, 4, seed=1),
            _balanced_dataset(4, 4, seed=2),
        ),
        (
            _balanced_dataset(6, 4, seed=3),
            _balanced_dataset(4, 4, seed=4),
        ),
    ]

    def builder(**params):
        return create_model_instance("xgboost", **params)

    opt = KFoldRandomHyperparamOpt(model_builder=builder, max_iter=1)
    metric = dc.metrics.Metric(dc.metrics.roc_auc_score)

    best_models, best_hyperparams, all_scores = opt.hyperparam_kfold_search(
        params_dict={"max_depth": [3], "n_tasks": [1]},
        folds=folds,
        metric=metric,
        use_max=True,
        logdir=str(tmp_path),
    )

    assert len(best_models) == len(folds)
    assert isinstance(best_hyperparams, dict)
    assert all_scores  # at least one hyperparameter combination scored
    assert all(0.0 <= score <= 1.0 for score in all_scores.values())
