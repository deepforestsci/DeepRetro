import os
import tempfile
import numpy as np
import deepchem as dc
from xgboost import XGBClassifier
from deepretro.featurizers.reactionstep import ReactionStepFeaturizer


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

FEATURIZER_MAP = {
    "molgraphconv": {
        "class": dc.feat.MolGraphConvFeaturizer,
        "type": "graph",
        "default_params": {}
    },
    "circularfingerprint": {
        "class": dc.feat.CircularFingerprint,
        "type": "1d_vector",
        "default_params": {"radius": 2, "size": 2048}
    },
    "reactionfeaturizer": {
        "class": ReactionStepFeaturizer,
        "type": "composite", 
        "default_params": {}
    }
}

def create_model_instance(model_type, **model_params):
    """
    Standalone utility function for reconstructing deepchem models anywhere.
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
