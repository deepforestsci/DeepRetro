"""XGBoost hallucination classifier built on DeepChem's GBDTModel.

Provides a single class that handles training, evaluation, threshold
optimisation, single-reaction prediction, and persistence using
DeepChem APIs end-to-end.  ``GBDTModel`` wraps an ``XGBClassifier``
and adds automatic early-stopping with an 80/20 internal split.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from deepchem.data import Dataset, NumpyDataset
from deepchem.metrics import Metric, accuracy_score, f1_score, roc_auc_score
from deepchem.models import GBDTModel
from rdkit import Chem
from xgboost import XGBClassifier

from deepretro.featurizers import ReactionStepFeaturizer
from deepretro.utils.metrics import find_optimal_threshold


# Helper functions (kept outside the class per DeepChem convention)

def probability_scores(dataset: Dataset, model) -> dict[str, float]:
    """Compute ROC-AUC and optimal threshold from probabilities.

    Parameters
    ----------
    dataset : Dataset
        Labelled dataset with ``y`` ground-truth.
    model : XGBClassifier
        Fitted sklearn-compatible model with ``predict_proba``.

    Returns
    -------
    dict
        Keys: ``roc_auc``, ``optimal_threshold``, ``optimal_f1``.
    """
    y_true = dataset.y.flatten()
    probabilities = model.predict_proba(dataset.X)[:, 1]
    auc = roc_auc_score(y_true, probabilities)
    opt_thr, opt_f1 = find_optimal_threshold(y_true, probabilities)
    return {
        "roc_auc": auc,
        "optimal_threshold": opt_thr,
        "optimal_f1": opt_f1,
    }


def predict_single_reaction(
    clf: "HallucinationClassifier",
    product_smiles: str,
    reactants_smiles: str,
) -> dict[str, Any]:
    """Predict whether a single reaction step is hallucinated.

    This is a module-level helper so it can be used independently of the
    class method.

    Parameters
    ----------
    clf : HallucinationClassifier
        A fitted classifier instance.
    product_smiles : str
        SMILES of the target product.
    reactants_smiles : str
        SMILES of the proposed reactants (dot-separated).

    Returns
    -------
    result : dict
        Keys: ``is_hallucination`` (bool), ``probability`` (float).
        On invalid SMILES an ``error`` key is added instead.

    Examples
    --------
    >>> from deepretro.models import HallucinationClassifier
    >>> clf = HallucinationClassifier()
    >>> clf.load("saved_model/")                              # doctest: +SKIP
    >>> predict_single_reaction(clf, "CCO", "CC.O")           # doctest: +SKIP
    {'is_hallucination': False, 'probability': 0.12}
    >>> predict_single_reaction(clf, "GARBAGE", "CC.O")       # doctest: +SKIP
    {'error': 'Invalid SMILES', 'is_hallucination': None, 'probability': None}
    """
    if Chem.MolFromSmiles(product_smiles) is None or \
       Chem.MolFromSmiles(reactants_smiles) is None:
        return {
            "error": "Invalid SMILES",
            "is_hallucination": None,
            "probability": None,
        }

    if clf.featurizer is None:
        clf.featurizer = ReactionStepFeaturizer()

    X = clf.featurizer.featurize([(product_smiles, reactants_smiles)])
    ds = NumpyDataset(X=X)
    probability = clf.predict_probability(ds)[0]
    return {
        "is_hallucination": bool(probability >= clf.threshold),
        "probability": float(probability),
    }


class HallucinationClassifier(GBDTModel):
    """
    Binary classifier for detecting hallucinated retrosynthesis reactions.

    Inherits from DeepChem's ``GBDTModel`` which wraps an
    ``XGBClassifier`` and adds automatic early-stopping via an
    internal 80/20 train/validation split.

    Training data
    -------------
    Prepare a CSV with at least these columns:

    * ``product`` — SMILES of the target product.
    * ``reactants`` — SMILES of proposed reactants (dot-separated for
      multiple reactants, e.g. ``"CC.O"``).
    * ``label`` — ``1`` if the reaction is hallucinated, ``0`` if real.

    Then load and train:

    .. code-block:: python

       from deepretro.data import ReactionDataLoader, stratified_split
       from deepretro.models import HallucinationClassifier

       loader = ReactionDataLoader()
       ds = loader.create_dataset("path/to/your_dataset.csv")
       train, valid, test = stratified_split(ds)

       clf = HallucinationClassifier(model_dir="my_models/")
       clf.fit(train)
       scores = clf.evaluate(test)
       print(scores)

    The trained model is persisted via DeepChem's standard joblib
    serialisation.  To reload later::

       clf = HallucinationClassifier(model_dir="my_models/")
       clf.load("my_models/")

    Parameters
    ----------
    model_dir : str, optional
        Directory for DeepChem model checkpoints.  If ``None``, a
        temporary directory is used (see ``deepchem.models.Model``).
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

    def __init__(self, model_dir: str | None = None,
                 early_stopping_rounds: int = 50,
                 **xgb_kwargs: Any) -> None:
        params = {**self._DEFAULT_XGB, **xgb_kwargs}
        xgb = XGBClassifier(**params)
        super().__init__(
            model=xgb,
            model_dir=model_dir,
            early_stopping_rounds=early_stopping_rounds,
            eval_metric="logloss",
        )
        self.threshold: float = 0.5
        self.featurizer: ReactionStepFeaturizer | None = None

    # Training

    def fit(self, train_dataset: Dataset) -> None:
        """
        Train the model on a DeepChem ``Dataset``.

        ``GBDTModel`` automatically performs an internal 80/20
        train/validation split for early stopping.  The model is
        auto-saved to ``model_dir`` after training.

        Parameters
        ----------
        train_dataset : Dataset
            Training data produced by ``deepretro.data.loader``.

        Examples
        --------
        >>> clf.fit(train_ds)  # doctest: +SKIP
        """
        super().fit(train_dataset)
        super().save()

    # Evaluation

    def evaluate(self, test_dataset: Dataset, metrics=None) -> dict[str, float]:
        """
        Evaluate using DeepChem ``Metric`` objects.

        Returns label-based metrics, plus probability-based ROC-AUC
        and the optimal threshold.  Updates ``self.threshold`` to the
        optimal value and auto-saves the model state.

        Parameters
        ----------
        test_dataset : Dataset
            Held-out test data.
        metrics : list of dc.metrics.Metric, optional
            Label-based metrics to compute.  If ``None``, defaults to:
            ``[Metric(accuracy_score, name="accuracy"),``
            ``Metric(f1_score, name="f1")]``.
            Any ``sklearn.metrics`` function that accepts ``(y_true, y_pred)``
            can be wrapped with ``dc.metrics.Metric``, e.g.
            ``Metric(precision_score, name="precision")``.

        Returns
        -------
        scores : dict
            Contains each requested metric name (or ``accuracy``/``f1`` when
            defaults are used), plus ``roc_auc``, ``optimal_threshold``,
            and ``optimal_f1``.

        Examples
        --------
        >>> scores = clf.evaluate(test_ds)  # doctest: +SKIP
        >>> scores["roc_auc"]               # doctest: +SKIP
        0.92
        """
        if metrics is None:
            metrics = [
                Metric(accuracy_score, name="accuracy"),
                Metric(f1_score, name="f1"),
            ]
        label_scores = super().evaluate(test_dataset, metrics)

        prob_scores = probability_scores(test_dataset, self.model)
        self.threshold = prob_scores["optimal_threshold"]

        # Auto-save with updated threshold
        self.save(self.model_dir)

        scores = dict(label_scores)
        scores.update(prob_scores)
        return scores

    # Prediction

    def predict_probability(self, dataset: Dataset) -> np.ndarray:
        """
        Return hallucination probabilities for each sample.

        Parameters
        ----------
        dataset : Dataset
            Data to score.

        Returns
        -------
        probabilities : np.ndarray, shape (n_samples,)
            Probability of the positive class (hallucination).
        """
        return self.model.predict_proba(dataset.X)[:, 1]

    def predict_with_threshold(self, dataset: Dataset) -> tuple[np.ndarray, np.ndarray]:
        """
        Predict binary labels using the current threshold.

        Unlike the inherited ``predict()`` (which returns raw model
        output), this applies ``self.threshold`` to produce binary labels.

        Parameters
        ----------
        dataset : Dataset
            Data to classify.

        Returns
        -------
        labels : np.ndarray, shape (n_samples,)
            Binary predictions (0 or 1).
        probabilities : np.ndarray, shape (n_samples,)
            Hallucination probabilities.
        """
        probabilities = self.predict_probability(dataset)
        return (probabilities >= self.threshold).astype(int), probabilities

    def predict_single(self, product_smiles: str, reactants_smiles: str) -> dict[str, Any]:
        """Thin wrapper around :func:`predict_single_reaction`."""
        return predict_single_reaction(self, product_smiles, reactants_smiles)

    # Persistence

    def save(self, save_dir: str) -> None:
        """
        Save model and threshold via DeepChem's joblib persistence.

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

        self.model.optimal_threshold_ = self.threshold
        self.model_dir = str(save_path)
        super().save()

    def load(self, save_dir: str) -> None:
        """
        Reload a previously saved model.

        Parameters
        ----------
        save_dir : str
            Directory containing saved artifacts.

        Examples
        --------
        >>> clf = HallucinationClassifier()
        >>> clf.load("saved_model/")  # doctest: +SKIP
        """
        self.model_dir = str(Path(save_dir))
        self.reload()
        self.threshold = getattr(self.model, "optimal_threshold_", 0.5)
        self.featurizer = ReactionStepFeaturizer()