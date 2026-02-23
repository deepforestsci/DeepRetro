"""XGBoost hallucination classifier using DeepChem's GBDTModel.

Provides a single class that handles training, evaluation, threshold
optimisation, single-reaction prediction, and persistence — using
DeepChem APIs end-to-end.  ``GBDTModel`` wraps an ``XGBClassifier``
and adds automatic early-stopping with an 80/20 internal split.
"""

import json
from pathlib import Path

import numpy as np
from deepchem.data import NumpyDataset
from deepchem.metrics import Metric, accuracy_score, f1_score, roc_auc_score
from deepchem.models import GBDTModel
from xgboost import XGBClassifier

from deepretro.featurizers import ReactionStepFeaturizer
from deepretro.utils.metrics import find_optimal_threshold


class HallucinationClassifier:
    """
    Binary classifier for detecting hallucinated retrosynthesis reactions.

    Wraps ``XGBClassifier`` inside DeepChem's ``GBDTModel`` which adds
    automatic early-stopping via an internal 80/20 train/validation
    split.  The training pipeline uses DeepChem datasets, splitters,
    and metrics end-to-end.

    Parameters
    ----------
    model_dir : str, optional
        Directory for DeepChem model checkpoints.  Default ``"model_dir"``.
    early_stopping_rounds : int, optional
        Rounds for early stopping during ``fit()``.  Default ``50``.
    **xgb_kwargs
        Forwarded to the underlying ``XGBClassifier`` via ``GBDTModel``.
        Defaults are tuned for the hallucination detection task.

    Examples
    --------
    >>> from deepretro.models import HallucinationClassifier
    >>> clf = HallucinationClassifier()
    >>> clf.threshold
    0.5
    """

    # Default XGBoost hyper-parameters
    _DEFAULT_XGB = dict(
        max_depth=6,
        learning_rate=0.05,
        n_estimators=300,
        subsample=0.8,
        colsample_bytree=0.8,
        min_child_weight=3,
        gamma=0.1,
        random_state=42,
    )

    def __init__(self, model_dir="model_dir", early_stopping_rounds=50,
                 **xgb_kwargs):
        params = {**self._DEFAULT_XGB, **xgb_kwargs}
        xgb = XGBClassifier(**params)
        self.dc_model = GBDTModel(
            model=xgb,
            model_dir=model_dir,
            early_stopping_rounds=early_stopping_rounds,
            eval_metric="logloss",
        )
        self.threshold: float = 0.5
        self.featurizer: ReactionStepFeaturizer | None = None

    # Training

    def fit(self, train_dataset):
        """
        Train the model on a DeepChem ``NumpyDataset``.

        ``GBDTModel`` automatically performs an internal 80/20
        train/validation split for early stopping.

        Parameters
        ----------
        train_dataset : NumpyDataset
            Training data produced by ``deepretro.data.loader``.

        Examples
        --------
        >>> clf.fit(train_ds)  # doctest: +SKIP
        """
        self.dc_model.fit(train_dataset)

    # Evaluation

    def evaluate(self, test_dataset):
        """
        Evaluate using DeepChem ``Metric`` objects.

        Returns label-based accuracy and F1, plus probability-based
        ROC-AUC and the optimal threshold.

        Parameters
        ----------
        test_dataset : NumpyDataset
            Held-out test data.

        Returns
        -------
        scores : dict
            Keys: ``roc_auc``, ``accuracy``, ``f1``,
            ``optimal_threshold``, ``optimal_f1``.

        Examples
        --------
        >>> scores = clf.evaluate(test_ds)  # doctest: +SKIP
        >>> scores["roc_auc"]               # doctest: +SKIP
        0.92
        """
        # Label-based metrics via DeepChem
        dc_metrics = [
            Metric(accuracy_score, name="accuracy"),
            Metric(f1_score, name="f1"),
        ]
        label_scores = self.dc_model.evaluate(test_dataset, dc_metrics)

        # Probability-based metrics (need predict_proba from underlying XGB)
        y_true = test_dataset.y.flatten()
        proba = self.dc_model.model.predict_proba(test_dataset.X)[:, 1]
        auc = roc_auc_score(y_true, proba)

        # Optimal threshold
        opt_thr, opt_f1 = find_optimal_threshold(y_true, proba)
        self.threshold = opt_thr

        return {
            "roc_auc": auc,
            "accuracy": label_scores["accuracy"],
            "f1": label_scores["f1"],
            "optimal_threshold": opt_thr,
            "optimal_f1": opt_f1,
        }

    # Prediction

    def predict_proba(self, dataset):
        """
        Return hallucination probabilities for each sample.

        Parameters
        ----------
        dataset : NumpyDataset
            Data to score.

        Returns
        -------
        proba : np.ndarray, shape (n_samples,)
            Probability of the positive class (hallucination).
        """
        return self.dc_model.model.predict_proba(dataset.X)[:, 1]

    def predict(self, dataset):
        """
        Predict binary labels using the current threshold.

        Parameters
        ----------
        dataset : NumpyDataset
            Data to classify.

        Returns
        -------
        labels : np.ndarray, shape (n_samples,)
            Binary predictions (0 or 1).
        proba : np.ndarray, shape (n_samples,)
            Hallucination probabilities.
        """
        proba = self.predict_proba(dataset)
        return (proba >= self.threshold).astype(int), proba

    def predict_single(self, product_smiles, reactants_smiles):
        """
        Predict whether a single reaction step is hallucinated.

        Parameters
        ----------
        product_smiles : str
            SMILES of the target product.
        reactants_smiles : str
            SMILES of the proposed reactants (dot-separated).

        Returns
        -------
        result : dict
            Keys: ``is_hallucination`` (bool), ``probability`` (float).

        Examples
        --------
        >>> from deepretro.models import HallucinationClassifier
        >>> clf = HallucinationClassifier()
        >>> clf.load("saved_model/")                     # doctest: +SKIP
        >>> clf.predict_single("CCO", "CC.O")            # doctest: +SKIP
        {'is_hallucination': False, 'probability': 0.12}
        """
        if self.featurizer is None:
            self.featurizer = ReactionStepFeaturizer()

        X = self.featurizer.featurize([(product_smiles, reactants_smiles)])
        ds = NumpyDataset(X=X)
        proba = self.predict_proba(ds)[0]
        return {
            "is_hallucination": bool(proba >= self.threshold),
            "probability": float(proba),
        }

    # Persistence

    def save(self, save_dir):
        """
        Save model, featurizer reference, and metadata.

        Parameters
        ----------
        save_dir : str
            Directory to write artifacts into.

        Examples
        --------
        >>> clf.save("saved_model/")  # doctest: +SKIP
        """
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)

        # DeepChem model checkpoint
        self.dc_model.model_dir = str(save_path)
        self.dc_model.save()

        # Metadata JSON
        meta = {
            "model_type": "XGBoost (via dc.models.GBDTModel)",
            "wrapper": "dc.models.GBDTModel",
            "optimal_threshold": self.threshold,
            "feature_dimension": self.featurizer.feature_dim if self.featurizer else None,
            "domain_features": True,
        }
        with open(save_path / "metadata.json", "w") as f:
            json.dump(meta, f, indent=2)

    def load(self, save_dir):
        """
        Reload a previously saved model and metadata.

        Parameters
        ----------
        save_dir : str
            Directory containing saved artifacts.

        Examples
        --------
        >>> clf = HallucinationClassifier()
        >>> clf.load("saved_model/")  # doctest: +SKIP
        """
        save_path = Path(save_dir)

        self.dc_model.model_dir = str(save_path)
        self.dc_model.reload()

        meta_file = save_path / "metadata.json"
        if meta_file.exists():
            with open(meta_file) as f:
                meta = json.load(f)
            self.threshold = meta.get("optimal_threshold", 0.5)

        self.featurizer = ReactionStepFeaturizer()