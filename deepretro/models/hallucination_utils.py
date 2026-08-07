"""Shared building blocks for the hallucination-detection ML pipeline.

This module holds the pieces that both the trainer and the checker rely on:

* :data:`MODEL_CONFIGS` -- the registry describing each supported model type
  (its DeepChem wrapper, base estimator, default hyperparameters, search
  space, and default featurizer).
* :func:`create_model_instance` -- a small factory that reconstructs a
  DeepChem model from a registry entry plus caller-supplied parameters.
* :class:`KFoldRandomHyperparamOpt` -- a DeepChem ``RandomHyperparamOpt``
  subclass that scores each hyperparameter combination by K-Fold
  cross-validation instead of a single train/validation split.
"""

import os
import tempfile
import numpy as np
import deepchem as dc
from xgboost import XGBClassifier


# Unified Infrastructure Model Configuration Mapping (goes with trainer)
MODEL_CONFIGS = {
    'xgboost': {
        'model_class': dc.models.SklearnModel,
        'base_estimator': XGBClassifier,
        'is_sklearn': True,
        'default_params': {'max_depth': 3},
        # 'param_space': {'max_depth': [3, 5, 7], 'learning_rate': [0.01, 0.05, 0.1]},
        'param_space' : {'eta': [0.01, 0.05, 0.1, 0.2, 0.3],
                        'max_depth': [3, 4, 5, 6, 7, 8, 10],
                        'gamma': [0.0, 0.1, 0.2, 0.5, 1.0, 5.0],
                        'min_child_weight': [1, 3, 5, 10],
                        'subsample': [0.6, 0.7, 0.8, 0.9, 1.0],
                        'colsample_bytree': [0.6, 0.7, 0.8, 0.9, 1.0],
                        'reg_lambda': [0.0, 0.1, 1.0, 5.0, 10.0],
                        'reg_alpha': [0.0, 0.1, 1.0, 5.0, 10.0]},
        'default_feat': 'reactionfeaturizer',
        'default_feat_params': {'base_featurizer': 'circularfingerprint', 'use_domain_features': True}
    },
#     'gcn': {
#         'model_class': dc.models.GCNModel,
#         'is_sklearn': False,
#         'default_params': {'batch_size': 32, 'learning_rate': 0.001},
#         'param_space': {'batch_size': [16, 32, 64], 'learning_rate': [1e-4, 1e-3, 1e-2]},
#         'default_feat': 'reactionfeaturizer',
#         'default_feat_params': {'base_featurizer': 'molgraphconv', 'use_domain_features': False}
#     }
}


def create_model_instance(model_type, **model_params):
    """Reconstruct a DeepChem model from the registry and caller parameters.

    This is a standalone factory so a model can be rebuilt anywhere (training,
    tuning, or inference) without threading a config object through the call
    stack. The infrastructure keys ``model_dir`` and ``n_tasks`` are consumed
    here; every remaining keyword is forwarded to the underlying estimator.

    Parameters
    ----------
    model_type : str
        Registry key selecting the model configuration (case-insensitive),
        e.g. ``"xgboost"``. Must be present in :data:`MODEL_CONFIGS`.
    **model_params
        Model hyperparameters. ``model_dir`` (str or None) and ``n_tasks``
        (int, default ``1``) are popped as infrastructure settings; all other
        keywords are passed to the base estimator (for sklearn models) or to
        the DeepChem model constructor.

    Returns
    -------
    deepchem.models.Model
        A ready-to-fit model instance. For ``"xgboost"`` this is a
        :class:`deepchem.models.SklearnModel` wrapping an ``XGBClassifier``.

    Raises
    ------
    KeyError
        If ``model_type`` is not registered in :data:`MODEL_CONFIGS`.

    Examples
    --------
    >>> from deepretro.models.hallucination_utils import create_model_instance
    >>> model = create_model_instance("xgboost", max_depth=3, n_tasks=1)
    >>> type(model).__name__
    'SklearnModel'
    >>> model.model.get_params()["max_depth"]
    3
    """
    model_dir = model_params.pop('model_dir', None)
    n_tasks = model_params.pop('n_tasks', 1)

    config = MODEL_CONFIGS[model_type.lower()]
    model_class = config['model_class']
    is_sklearn = config.get('is_sklearn', False)

    if is_sklearn:
        base_estimator_class = config['base_estimator']
        base_estimator = base_estimator_class(**model_params)
        return model_class(base_estimator, model_dir=model_dir)
    else:
        kwargs = model_params.copy()
        kwargs['n_tasks'] = n_tasks
        kwargs['mode'] = 'classification'
        kwargs['model_dir'] = model_dir
        return model_class(**kwargs)


class KFoldRandomHyperparamOpt(dc.hyper.RandomHyperparamOpt):
    """
    Child class of DeepChem's RandomHyperparamOpt that introduces K-Fold cross-validation support.
    """

    def hyperparam_kfold_search(self, params_dict, folds, metric, use_max=True, logdir=None, **kwargs):
        """Search hyperparameters, scoring each combination by K-Fold CV.

        For every randomly sampled hyperparameter combination, a model is
        trained and evaluated on each fold and the fold scores are averaged.
        The combination with the best average validation score wins.

        Parameters
        ----------
        params_dict : dict
            Mapping of hyperparameter name to the list of candidate values to
            sample from.
        folds : Iterable[tuple[deepchem.data.Dataset, deepchem.data.Dataset]]
            Sequence of ``(train_dataset, valid_dataset)`` pairs, one per fold.
        metric : deepchem.metrics.Metric
            Metric used to score each fold's validation set.
        use_max : bool, optional
            If ``True`` (default) higher metric values are better; if ``False``
            lower values are better.
        logdir : str, optional
            Directory to persist per-fold model artifacts. When ``None`` a
            temporary directory is used for each fold.
        **kwargs
            Extra keyword arguments forwarded to each model's ``fit`` call
            (silently dropped if the model's ``fit`` does not accept them).

        Returns
        -------
        tuple[list, dict, dict]
            ``(best_models, best_hyperparams, all_scores)`` where
            ``best_models`` holds one trained model per fold for the winning
            combination, ``best_hyperparams`` is that combination, and
            ``all_scores`` maps each combination's filename-encoded key to its
            average validation score.
        """
        hyperparameter_combs = self.generate_random_hyperparam_values(params_dict, self.max_iter)

        if use_max:
            best_validation_score = -np.inf
        else:
            best_validation_score = np.inf

        best_hyperparams = None
        best_models = None
        all_scores = {}

        for ind, model_params in enumerate(hyperparameter_combs):
            hp_str = dc.hyper.base_classes._convert_hyperparam_dict_to_filename(model_params)
            fold_scores = []
            current_fold_models = []

            for fold_idx, (train_dataset, valid_dataset) in enumerate(folds):
                if logdir is not None:
                    model_dir = os.path.join(logdir, hp_str, f"fold_{fold_idx}")
                    os.makedirs(model_dir, exist_ok=True)
                else:
                    model_dir = tempfile.mkdtemp()

                model_params['model_dir'] = model_dir
                model = self.model_builder(**model_params)

                try:
                    model.fit(train_dataset, **kwargs)
                except TypeError:
                    model.fit(train_dataset)

                current_fold_models.append(model)
                multitask_scores = model.evaluate(valid_dataset, [metric])
                fold_scores.append(multitask_scores[metric.name])

            avg_valid_score = np.mean(fold_scores)
            all_scores[hp_str] = avg_valid_score

            if (use_max and avg_valid_score >= best_validation_score) or (not use_max and avg_valid_score <= best_validation_score):
                best_validation_score = avg_valid_score
                best_hyperparams = model_params
                best_models = current_fold_models

        return best_models, best_hyperparams, all_scores
