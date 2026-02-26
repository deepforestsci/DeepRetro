"""Dataset loading pipeline using DeepChem data structures.

Converts a reaction CSV (product, reactants, label) into DeepChem
``NumpyDataset`` objects with stratified train/valid/test splits.
"""

import warnings

import numpy as np
import pandas as pd
from deepchem.data import NumpyDataset
from deepchem.splits import SingletaskStratifiedSplitter

from deepretro.featurizers import ReactionStepFeaturizer


def load_reaction_csv(
    csv_path,
    product_col="product",
    reactants_col="reactants",
    label_col="label",
):
    """
    Read a reaction CSV into parallel lists.

    Parameters
    ----------
    csv_path : str
        Path to CSV file.
    product_col : str, optional
        Column name for product SMILES. Default ``"product"``.
    reactants_col : str, optional
        Column name for reactant SMILES. Default ``"reactants"``.
    label_col : str, optional
        Column name for binary labels. Default ``"label"``.

    Returns
    -------
    products : list of str
    reactants : list of str
    labels : list of int

    Examples
    --------
    >>> from deepretro.data import load_reaction_csv
    >>> products, reactants, labels = load_reaction_csv("data/dataset.csv")
    """
    df = pd.read_csv(csv_path)
    return (
        df[product_col].tolist(),
        df[reactants_col].tolist(),
        df[label_col].tolist(),
    )


def featurize_reactions(products, reactants, featurizer=None):
    """
    Featurize product-reactant pairs using ``ReactionStepFeaturizer``.

    Parameters
    ----------
    products : list of str
        Product SMILES strings.
    reactants : list of str
        Reactant SMILES strings (dot-separated when multiple).
    featurizer : ReactionStepFeaturizer, optional
        Pre-configured featurizer.  A default one (radius=2, size=2048,
        domain features on) is created when ``None``.

    Returns
    -------
    X : np.ndarray, shape (n_samples, feature_dim)
        Feature matrix.
    featurizer : ReactionStepFeaturizer
        The featurizer instance (useful for saving alongside the model).

    Examples
    --------
    >>> from deepretro.data import featurize_reactions
    >>> X, feat = featurize_reactions(["CCO"], ["CC.O"])
    >>> X.shape[1]
    4111
    """
    if featurizer is None:
        featurizer = ReactionStepFeaturizer()
    pairs = list(zip(products, reactants))
    X = featurizer.featurize(pairs)
    return X, featurizer


def create_dataset(X, y):
    """
    Wrap a feature matrix and labels into a DeepChem ``NumpyDataset``.

    Parameters
    ----------
    X : np.ndarray, shape (n_samples, n_features)
        Feature matrix.
    y : array-like, shape (n_samples,)
        Binary labels.

    Returns
    -------
    dataset : NumpyDataset

    Examples
    --------
    >>> import numpy as np
    >>> from deepretro.data import create_dataset
    >>> ds = create_dataset(np.zeros((5, 10)), [0, 1, 0, 1, 0])
    >>> len(ds)
    5
    """
    y = np.array(y).reshape(-1, 1)
    nan_mask = np.isnan(X).all(axis=1)
    if nan_mask.any():
        warnings.warn(
            f"Dropped {nan_mask.sum()} rows with NaN features (invalid SMILES)."
        )
        X = X[~nan_mask]
        y = y[~nan_mask]
    return NumpyDataset(X=X, y=y)


def split_dataset(dataset, frac_train=0.7, frac_valid=0.15, frac_test=0.15, seed=42):
    """
    Stratified split into train / valid / test sets.

    Uses DeepChem's ``SingletaskStratifiedSplitter`` to preserve class
    balance across all three splits.

    Parameters
    ----------
    dataset : NumpyDataset
        Full dataset to split.
    frac_train : float, optional
        Training fraction. Default 0.7.
    frac_valid : float, optional
        Validation fraction. Default 0.15.
    frac_test : float, optional
        Test fraction. Default 0.15.
    seed : int, optional
        Random seed. Default 42.

    Returns
    -------
    train_ds : NumpyDataset
    valid_ds : NumpyDataset
    test_ds : NumpyDataset

    Examples
    --------
    >>> import numpy as np
    >>> from deepretro.data import create_dataset, split_dataset
    >>> ds = create_dataset(np.random.rand(100, 10), [0]*50 + [1]*50)
    >>> train, valid, test = split_dataset(ds)
    >>> len(train) + len(valid) + len(test) == 100
    True
    """
    splitter = SingletaskStratifiedSplitter()
    train_inds, valid_inds, test_inds = splitter.split(
        dataset,
        frac_train=frac_train,
        frac_valid=frac_valid,
        frac_test=frac_test,
        seed=seed,
    )
    return (
        dataset.select(train_inds.astype(int)),
        dataset.select(valid_inds.astype(int)),
        dataset.select(test_inds.astype(int)),
    )

