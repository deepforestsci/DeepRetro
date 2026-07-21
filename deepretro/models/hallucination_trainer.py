import os
import json
import shutil
import numpy as np
import deepchem as dc
from deepretro.featurizers.reactionstep import FEATURIZER_MAP
from deepretro.models.hallucination_utils import (
    MODEL_CONFIGS,
    create_model_instance,
    KFoldRandomHyperparamOpt,
)
from sklearn.metrics import (
    precision_recall_curve,
    accuracy_score,
    precision_score,
    recall_score,
)
from deepretro.data.loader import ReactionDataLoader
from typing import Any, Optional, Dict, Tuple, List


class HallucinationTrainer:
    """
    A DeepChem-based ML infrastructure for training models to detect hallucinations.
    """

    def __init__(
        self,
        trainer_dir: str,
        n_tasks: int = 1,
        model_type: Optional[str] = None,
        model_path: Optional[str] = None,
    ) -> None:
        """
        Initializes a hallucination trainer.

        Parameters
        ----------
        trainer_dir : str
            Directory used for storing model artifacts.
        n_tasks : int
            Number of prediction tasks, default=1
        model_type : str, optional
            Name of the model configuration to initialize.
        model_path : str, optional
            Path to a previously saved model configuration.
        """
        self.trainer_dir = trainer_dir
        os.makedirs(self.trainer_dir, exist_ok=True)

        # Initialize attributes
        self.model_type = None
        self.model_params = {}
        self.feat_name = None
        self.threshold = 0.5
        self.param_space = {}
        self.model_dir = None
        self.model = None
        self.n_tasks = n_tasks
        self.featurizer = None
        self.feat_params = None

        if model_path is not None:
            config_path = os.path.join(model_path, "config.json")
            if not os.path.exists(config_path):
                raise FileNotFoundError(f"No config.json found in {model_path}")

            with open(config_path, "r") as f:
                config = json.load(f)

            self.model_type = config.get("model_type")
            self.model_params = config.get("model_params", {})
            self.feat_name = config.get("feat_name")
            self.feat_params = config.get("feat_params", {})
            self.threshold = config.get("threshold", 0.5)
            self.param_space = config.get("param_space", {})
            self.n_tasks = config.get("n_tasks", self.n_tasks)

            self.model_dir = os.path.join(self.trainer_dir, f"{self.model_type}_model")
            os.makedirs(self.model_dir, exist_ok=True)

            if os.path.abspath(model_path) != os.path.abspath(self.model_dir):
                for item in os.listdir(model_path):
                    s = os.path.join(model_path, item)
                    d = os.path.join(self.model_dir, item)
                    if os.path.isdir(s):
                        shutil.copytree(s, d, dirs_exist_ok=True)
                    else:
                        shutil.copy2(s, d)

        elif model_type is not None:
            self.model_type = model_type.lower()
            if self.model_type not in MODEL_CONFIGS:
                raise ValueError(
                    f"Unknown model_type: {self.model_type}. Available: {list(MODEL_CONFIGS.keys())}"
                )

            config = MODEL_CONFIGS[self.model_type]
            self.model_params = config["default_params"]
            self.param_space = config["param_space"]
            self.feat_name = config["default_feat"]

            # Adopt model default featurizer parameter mappings if none are supplied externally
            if self.feat_params is None:
                self.feat_params = config["default_feat_params"].copy()

            self.model_dir = os.path.join(self.trainer_dir, f"{self.model_type}_model")
            os.makedirs(self.model_dir, exist_ok=True)
        else:
            raise ValueError(
                "You must provide either a 'model_type' (for a new model) or a 'model_path' (to load an existing one)."
            )

        # Dynamically build ReactionStepFeaturizer wrapper using the model-specific base configuration mappings
        feat_config = FEATURIZER_MAP[self.feat_name.lower()]
        self.featurizer = feat_config["class"](**self.feat_params)

    def _build_model(self, n_tasks: int) -> None:
        """
        Instantiates the configured model and restores saved weights when available.

        Parameters
        ----------
        n_tasks : int
            Number of prediction tasks.
        """
        params = self.model_params.copy()
        params["model_dir"] = self.model_dir
        params["n_tasks"] = n_tasks

        self.model = create_model_instance(self.model_type, **params)

        if os.path.exists(self.model_dir) and os.listdir(self.model_dir):
            try:
                self.model.restore()
            except Exception:
                self.model.reload()

    def load_dataset(
        self,
        train_csv: str,
        test_csv: str,
        product_col: str = "product",
        reactants_col: str = "reactants",
        label_col: str = "label",
    ) -> Tuple[dc.data.Dataset, dc.data.Dataset]:
        """
        Loads and featurizes train and test datasets.

        Parameters
        ----------
        train_csv: str
            Path to the train CSV file.
        test_csv: str
            Path to the test CSV file.
        product_col: str
            Product SMILES column name, default="product"
        reactants_col: str
            Reactants SMILES column name, default="reactants"
        label_col: str
            Label column name, default="label"

        Returns
        -------
        tuple[Dataset, Dataset]
            Train and test datasets.
        """
        if self.n_tasks != 1:
            raise ValueError(
                f"Configured n_tasks ({self.n_tasks}) does not match the number of label columns ({1})."
            )

        loader = ReactionDataLoader(
            featurizer=self.featurizer,
            product_col=product_col,
            reactants_col=reactants_col,
            label_col=label_col,
        )

        print(f"Loading training data from {train_csv}")
        train_dataset = loader.create_dataset(train_csv)

        print(f"Loading testing data from {test_csv}")
        test_dataset = loader.create_dataset(test_csv)

        return train_dataset, test_dataset

    def optimize_threshold(
        self,
        model: dc.models.Model,
        valid_dataset: dc.data.Dataset,
    ) -> float:
        """
        Computes the F1 scores across precision-recall thresholds and returns the best threshold.

        Parameters
        ----------
        model: dc.models.Model
            Trained DeepChem model
        valid_dataset: dc.data.Dataset
            Validation dataset

        Returns
        -------
        float
            Threshold that maximizes F1 score.
        """
        y_true = valid_dataset.y.flatten()
        y_pred = model.predict(valid_dataset)

        # Handles both (N, 1, 2) and (N, 2) shapes cleanly without corruption
        if len(y_pred.shape) == 3 and y_pred.shape[2] == 2:
            y_pred_proba = y_pred[:, :, 1].flatten()
        elif len(y_pred.shape) == 2 and y_pred.shape[1] == 2:
            y_pred_proba = y_pred[:, 1].flatten()
        else:
            y_pred_proba = y_pred.flatten()

        precisions, recalls, thresholds = precision_recall_curve(y_true, y_pred_proba)

        f1_scores = np.divide(
            2 * (precisions * recalls),
            (precisions + recalls),
            out=np.zeros_like(precisions),
            where=(precisions + recalls) != 0,
        )

        best_idx = np.argmax(f1_scores)

        if best_idx == len(thresholds):
            best_threshold = thresholds[-1]
        else:
            best_threshold = thresholds[best_idx]

        return best_threshold

    def parameter_tuning(
        self,
        train_dataset: dc.data.Dataset,
        n_trials: int = 10,
        k_folds: int = 1,
    ) -> tuple[dict[str, Any], float]:
        """
        Performs hyperparameter optimization.

        Parameters
        ----------
        train_dataset : dc.data.Dataset
            Dataset used for tuning.
        n_trials : int
            Number of hyperparameter trials, default=10
        k_folds : int
            Number of cross-validation folds, default=1

        Returns
        -------
        tuple[dict[str, Any], float]
            Best parameter set and optimal threshold.
        """

        if not self.param_space:
            print("No parameter space defined. Skipping tuning.")
            return self.model_params, self.threshold

        tune_space = self.param_space.copy()
        tune_space["n_tasks"] = [self.n_tasks]

        eval_metric = dc.metrics.Metric(dc.metrics.prc_auc_score)

        def lambda_model_builder(**params):
            return create_model_instance(self.model_type, **params)

        if k_folds <= 1:
            print(
                f"Starting single-split random hyperparameter tuning with {n_trials} trials..."
            )
            splitter = dc.splits.RandomStratifiedSplitter()
            tune_train, tune_valid = splitter.train_test_split(
                train_dataset, frac_train=0.8
            )

            optimizer = dc.hyper.RandomHyperparamOpt(
                model_builder=lambda_model_builder, max_iter=n_trials
            )

            best_model, best_params, all_scores = optimizer.hyperparam_search(
                params_dict=tune_space,
                train_dataset=tune_train,
                valid_dataset=tune_valid,
                metric=eval_metric,
                use_max=True,
            )

            best_threshold = self.optimize_threshold(best_model, tune_valid)

        else:
            print(
                f"Starting K-Fold ({k_folds}) random hyperparameter tuning with {n_trials} trials..."
            )
            splitter = dc.splits.RandomStratifiedSplitter()
            folds = list(splitter.k_fold_split(train_dataset, k_folds))

            optimizer = KFoldRandomHyperparamOpt(
                model_builder=lambda_model_builder, max_iter=n_trials
            )

            best_models, best_params, all_scores = optimizer.hyperparam_kfold_search(
                params_dict=tune_space, folds=folds, metric=eval_metric, use_max=True
            )
            print(
                "Calculating average optimal threshold across folds for best parameters..."
            )
            fold_thresholds = []
            for fold_model, (_, valid_dataset) in zip(best_models, folds):
                thresh = self.optimize_threshold(fold_model, valid_dataset)
                fold_thresholds.append(thresh)

            best_threshold = np.mean(fold_thresholds)

        best_score = all_scores.get(
            dc.hyper.base_classes._convert_hyperparam_dict_to_filename(best_params), 0.0
        )
        best_params.pop("n_tasks", None)
        best_params.pop("model_dir", None)

        print("\n--- Tuning Complete ---")
        print(f"Best Params: {best_params}")
        print(f"Best Threshold: {best_threshold:.4f} (PRC-AUC Score: {best_score:.4f})")

        return best_params, best_threshold

    def evaluate(
        self,
        dataset: dc.data.Dataset,
        metrics: List[dc.metrics.Metric] | None = None,
    ) -> Dict[str, float]:
        """
        Evaluates the model on a given dataset using base metrics and threshold-based metrics.

        Parameters
        ----------
        dataset : dc.data.Dataset
            Dataset to evaluate.
        metrics : list[Metric]
            DeepChem metrics to compute, Optional

        Returns
        -------
        dict[str, float]
            Evaluation metrics and threshold-based statistics.
        """

        if self.model is None:
            raise RuntimeError(
                "Model has not been instantiated. Ensure _build_model is triggered before evaluating."
            )

        if metrics is None:
            metrics = [
                dc.metrics.Metric(dc.metrics.roc_auc_score, name="ROC_AUC"),
                dc.metrics.Metric(dc.metrics.prc_auc_score, name="PRC_AUC"),
            ]

        print("Evaluating dataset...")
        scores = self.model.evaluate(dataset, metrics)

        y_true = dataset.y.flatten()
        y_pred = self.model.predict(dataset)

        if len(y_pred.shape) == 3 and y_pred.shape[2] == 2:
            y_pred_proba = y_pred[:, :, 1].flatten()
        elif len(y_pred.shape) == 2 and y_pred.shape[1] == 2:
            y_pred_proba = y_pred[:, 1].flatten()
        else:
            y_pred_proba = y_pred.flatten()

        y_pred_bin = (y_pred_proba >= self.threshold).astype(int)

        scores["threshold_accuracy"] = accuracy_score(y_true, y_pred_bin)
        scores["threshold_precision"] = precision_score(
            y_true, y_pred_bin, zero_division=0
        )
        scores["threshold_recall"] = recall_score(y_true, y_pred_bin, zero_division=0)

        print(f"Evaluation Scores: {scores}")
        return scores

    def save_config(self):
        """
        Saves the current trainer configuration to a JSON file in the model directory.
        """
        os.makedirs(self.model_dir, exist_ok=True)
        config_path = os.path.join(self.model_dir, "config.json")

        config = {
            "model_type": self.model_type,
            "model_params": self.model_params,
            "feat_name": self.feat_name,
            "feat_params": self.feat_params,
            "threshold": float(self.threshold),
            "param_space": self.param_space,
            "n_tasks": self.n_tasks,
        }

        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)

        print(f"Trainer configuration successfully saved to {config_path}")

    def train_model(
        self,
        train_dataset: dc.data.Dataset,
        test_dataset: dc.data.Dataset,
        tune_params: bool = False,
        n_trials: int = 10,
        k_folds: int = 1,
        **train_params: Any,
    ) -> tuple[dc.models.Model, dict[str, float]]:
        """
        Trains the final model on the provided dataset. Can optionally tune parameters first.

        Parameters
        ----------
        train_dataset: Dataset
            Train dataset
        test_dataset: Dataset
            Test dataset
        tune_params: bool
            To initiate Pre-Training parameter tuning, default=False
        n_trials : int
            Number of tuning trials, default=10
        k_folds: int
            Number of cross-validation folds, default=1
        **train_params: Any
            Additional model training parameters.

        Returns
        -------
        tuple[Model, dict[str, float]]
            Trained model and evaluation scores.

        """
        if tune_params:
            print("\n--- Initiating Pre-Training Parameter Tuning ---")
            best_params, best_threshold = self.parameter_tuning(
                train_dataset, n_trials=n_trials, k_folds=k_folds
            )
            self.model_params = best_params
            self.threshold = best_threshold

        print(f"\n--- Building and Training Final {self.model_type.upper()} Model ---")
        self._build_model(self.n_tasks)

        self.model.fit(train_dataset, **train_params)
        print("\n--- Final Training Complete. Evaluating on Test Dataset ---")
        test_scores = self.evaluate(test_dataset)

        print("\n--- Saving Model and Configuration ---")
        try:
            self.model.save()
        except NotImplementedError:
            raise NotImplementedError("Model save method not implemented.")

        self.save_config()
        return self.model, test_scores
