"""Dataset loading pipeline using DeepChem data structures.

Provides :class:`ReactionDataLoader`, a subclass of DeepChem's
:class:`~deepchem.data.DataLoader` that converts a reaction CSV
(product, reactants, label) into a ``DiskDataset`` with automatic
sharding, plus a convenience ``stratified_split`` function for
stratified train/valid/test splitting.
"""

from __future__ import annotations

import warnings
from typing import Iterator, List, Optional

import numpy as np
import pandas as pd
from deepchem.data import Dataset, DiskDataset
from deepchem.data.data_loader import DataLoader
from deepchem.splits import SingletaskStratifiedSplitter

from deepretro.featurizers import ReactionStepFeaturizer


class ReactionDataLoader(DataLoader):
    """Load a labelled reaction CSV into a DeepChem ``DiskDataset``.

    Inherits from :class:`~deepchem.data.DataLoader` and overrides
    :meth:`_get_shards` and :meth:`_featurize_shard` to handle paired
    SMILES columns (product + reactants).  The parent's
    ``create_dataset`` method is overridden to add automatic failed-row
    dropping (invalid SMILES) with a summary warning.

    Parameters
    ----------
    featurizer : ReactionStepFeaturizer or None, optional
        Pre-configured featurizer.  A default one (radius=2, size=2048,
        domain features off) is created when ``None``.
    product_col : str, optional
        Column name for product SMILES.  Default ``"product"``.
    reactants_col : str, optional
        Column name for reactant SMILES.  Default ``"reactants"``.
    label_col : str, optional
        Column name for binary labels.  Default ``"label"``.
    id_field : str or None, optional
        Column name to use as sample identifiers.  When ``None``
        (default), sequential integer IDs are generated per shard.
    log_every_n : int, optional
        DeepChem progress logging interval.  Default ``1000``.

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
        featurizer: ReactionStepFeaturizer | None = None,
        product_col: str = "product",
        reactants_col: str = "reactants",
        label_col: str = "label",
        id_field: str | None = None,
        log_every_n: int = 1000,
    ) -> None:
        """Store CSV column names and configure the DeepChem loader."""
        feat = featurizer or ReactionStepFeaturizer()
        super().__init__(
            tasks=[label_col],
            featurizer=feat,
            id_field=id_field,
            log_every_n=log_every_n,
        )
        self.product_col = product_col
        self.reactants_col = reactants_col
        self.label_col = label_col

    # DataLoader interface overrides

    def _get_shards(
        self,
        inputs: List[str],
        shard_size: Optional[int],
    ) -> Iterator[pd.DataFrame]:
        """Yield DataFrame chunks from one or more CSV files.

        Parameters
        ----------
        inputs : list of str
            Paths to CSV files.
        shard_size : int or None
            Rows per shard.  ``None`` reads the whole file at once.
        """
        for input_file in inputs:
            if shard_size is not None:
                for chunk in pd.read_csv(input_file, chunksize=shard_size):
                    yield chunk
            else:
                yield pd.read_csv(input_file)

    def _featurize_shard(
        self,
        shard: pd.DataFrame,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Featurize a shard of reaction pairs.

        Parameters
        ----------
        shard : pd.DataFrame
            A chunk of the input CSV.

        Returns
        -------
        X : np.ndarray
            Feature matrix (only valid rows).
        valid_inds : np.ndarray
            Integer indices of rows that were successfully featurized
            (all-NaN or empty feature rows are excluded).
        """
        products = shard[self.product_col].tolist()
        reactants = shard[self.reactants_col].tolist()
        pairs = list(zip(products, reactants))
        X = self.featurizer.featurize(pairs)

        # DeepChem represents failed graph featurization with an empty array.
        # Vector failures use an all-NaN row to retain the feature dimension.
        valid_inds = np.asarray(
            [
                index
                for index, row in enumerate(X)
                if not isinstance(row, np.ndarray)
                or (row.size > 0 and not np.isnan(row).all())
            ],
            dtype=int,
        )
        X = X[valid_inds]
        return X, valid_inds

    # create_dataset : customised to handle missing id_field & NaN warn

    # Unlike the generic loader, this implementation accepts CSV paths only.
    def create_dataset(  # ty: ignore[invalid-method-override]
        self,
        inputs: str | List[str],
        data_dir: Optional[str] = None,
        shard_size: Optional[int] = 1000,
    ) -> DiskDataset:
        """Read, featurize, and write a reaction CSV to a ``DiskDataset``.

        Follows the same contract as
        :meth:`~deepchem.data.DataLoader.create_dataset` but adds
        automatic dropping of rows where featurization fails (invalid
        SMILES), with a summary warning at the end.

        Parameters
        ----------
        inputs : str or list of str
            Path(s) to CSV file(s).
        data_dir : str or None, optional
            Directory to write shards into.  When ``None`` a temporary
            directory is used. Remove that directory when it is no longer needed.
        shard_size : int or None, optional
            Number of rows per shard.  Default ``1000``.

        Returns
        -------
        dataset : DiskDataset
            Disk-backed DeepChem dataset ready for splitting / training.

        Examples
        --------
        >>> loader = ReactionDataLoader()
        >>> ds = loader.create_dataset("data/dataset.csv")  # doctest: +SKIP
        """
        if not isinstance(inputs, list):
            inputs = [inputs]

        total_dropped = 0

        # Stream aligned features, labels, weights, and identifiers to disk.
        def shard_generator() -> Iterator[
            tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        ]:
            """Yield valid features, labels, weights, and IDs for each CSV shard."""
            nonlocal total_dropped
            for shard in self._get_shards(inputs, shard_size):
                X, valid_inds = self._featurize_shard(shard)

                dropped = len(shard) - len(valid_inds)
                if dropped > 0:
                    total_dropped += dropped

                # A custom vector width may still be unknown in an all-failed
                # shard. Do not persist empty shards with incompatible widths.
                if not len(valid_inds):
                    continue

                # Extract labels from task columns
                y = shard[self.label_col].values.reshape(-1, 1)
                y = y[valid_inds]

                w = np.ones_like(y, dtype=np.float32)

                # Use id_field if provided, else sequential IDs
                if self.id_field is not None and self.id_field in shard.columns:
                    ids = shard[self.id_field].values[valid_inds]
                else:
                    ids = np.arange(len(valid_inds))

                yield X, y, w, ids

        dataset = DiskDataset.create_dataset(shard_generator(), data_dir, self.tasks)

        if total_dropped > 0:
            warnings.warn(
                f"Dropped {total_dropped} rows with failed features (invalid SMILES).",
                stacklevel=2,
            )

        return dataset


def stratified_split(
    dataset: Dataset,
    frac_train: float = 0.7,
    frac_valid: float = 0.15,
    frac_test: float = 0.15,
    seed: int = 42,
) -> tuple[Dataset, Dataset, Dataset]:
    """Stratified split into train / valid / test sets.

    Uses DeepChem's ``SingletaskStratifiedSplitter`` to preserve class
    balance across all three splits.  Works with any DeepChem ``Dataset``
    subclass (``DiskDataset``, ``NumpyDataset``, etc.).

    Parameters
    ----------
    dataset : Dataset
        Full dataset to split (``DiskDataset`` or ``NumpyDataset``).
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
    train_ds : Dataset
    valid_ds : Dataset
    test_ds : Dataset

    Examples
    --------
    >>> import numpy as np
    >>> from deepchem.data import NumpyDataset
    >>> from deepretro.data import stratified_split
    >>> ds = NumpyDataset(X=np.random.rand(100, 10), y=np.array([0]*50 + [1]*50).reshape(-1,1))
    >>> train, valid, test = stratified_split(ds)
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
