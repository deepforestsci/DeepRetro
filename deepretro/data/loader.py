"""Dataset loading pipeline using DeepChem data structures.

Provides :class:`ReactionDataLoader`, a class-based loader modelled after
DeepChem's ``CSVLoader`` that converts a reaction CSV (product, reactants,
label) into ``NumpyDataset`` objects with optional stratified splitting.
"""

import warnings

import numpy as np
import pandas as pd
from deepchem.data import NumpyDataset
from deepchem.splits import SingletaskStratifiedSplitter

from deepretro.featurizers import ReactionStepFeaturizer


class ReactionDataLoader:
    """Load a two-column reaction CSV into a DeepChem ``NumpyDataset``.

    Mirrors the ``CSVLoader`` pattern: configuration is stored in
    ``__init__`` and the actual work happens in :meth:`create_dataset`.
    Unlike ``CSVLoader``, this loader accepts **two** SMILES columns
    (product + reactants) and returns an in-memory ``NumpyDataset``
    instead of a ``DiskDataset``.

    Parameters
    ----------
    featurizer : ReactionStepFeaturizer, optional
        Pre-configured featurizer.  A default one (radius=2, size=2048,
        domain features on) is created when ``None``.
    product_col : str, optional
        Column name for product SMILES.  Default ``"product"``.
    reactants_col : str, optional
        Column name for reactant SMILES.  Default ``"reactants"``.
    label_col : str, optional
        Column name for binary labels.  Default ``"label"``.

    Examples
    --------
    >>> from deepretro.data import ReactionDataLoader
    >>> loader = ReactionDataLoader()
    >>> ds = loader.create_dataset("data/dataset.csv")  # doctest: +SKIP
    >>> len(ds)                                          # doctest: +SKIP
    808
    """

    def __init__(
        self,
        featurizer=None,
        product_col="product",
        reactants_col="reactants",
        label_col="label",
    ):
        self.featurizer = featurizer or ReactionStepFeaturizer()
        self.product_col = product_col
        self.reactants_col = reactants_col
        self.label_col = label_col

    def create_dataset(self, csv_path):
        """Read, featurize, and wrap a reaction CSV into a ``NumpyDataset``.

        Rows where featurization fails (invalid SMILES) are automatically
        dropped with a warning.

        Parameters
        ----------
        csv_path : str
            Path to the CSV file.

        Returns
        -------
        dataset : NumpyDataset
            In-memory DeepChem dataset ready for splitting / training.

        Examples
        --------
        >>> loader = ReactionDataLoader()
        >>> ds = loader.create_dataset("data/dataset.csv")  # doctest: +SKIP
        """
        # 1. Read CSV
        df = pd.read_csv(csv_path)
        products = df[self.product_col].tolist()
        reactants = df[self.reactants_col].tolist()
        y = np.array(df[self.label_col]).reshape(-1, 1)

        # 2. Featurize
        pairs = list(zip(products, reactants))
        X = self.featurizer.featurize(pairs)

        # 3. Drop NaN rows (invalid SMILES)
        nan_mask = np.isnan(X).all(axis=1)
        if nan_mask.any():
            warnings.warn(
                f"Dropped {nan_mask.sum()} rows with NaN features "
                f"(invalid SMILES)."
            )
            X = X[~nan_mask]
            y = y[~nan_mask]

        return NumpyDataset(X=X, y=y)


def split_dataset(dataset, frac_train=0.7, frac_valid=0.15, frac_test=0.15, seed=42):
    """Stratified split into train / valid / test sets.

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
    >>> from deepchem.data import NumpyDataset
    >>> from deepretro.data import split_dataset
    >>> ds = NumpyDataset(X=np.random.rand(100, 10), y=np.array([0]*50 + [1]*50).reshape(-1,1))
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
