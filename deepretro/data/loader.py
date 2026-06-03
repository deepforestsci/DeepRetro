"""Dataset loading pipeline using DeepChem data structures.

Provides :class:`ReactionDataLoader`, a subclass of DeepChem's
:class:`~deepchem.data.DataLoader` that converts a reaction CSV
(product, reactants, label) into a ``DiskDataset`` with automatic
sharding, plus a convenience ``stratified_split`` function for
stratified train/valid/test splitting.
"""

from __future__ import annotations

import warnings
from typing import Any, Iterator, List, Optional

import numpy as np
import pandas as pd

_DEEPCHEM_IMPORT_ERROR: ModuleNotFoundError | None = None

try:
    from deepchem.data import Dataset, DiskDataset
    from deepchem.data.data_loader import DataLoader
    from deepchem.splits import SingletaskStratifiedSplitter
except ModuleNotFoundError as exc:
    _DEEPCHEM_IMPORT_ERROR = exc
    Dataset = Any  # type: ignore[assignment]
    DiskDataset = Any  # type: ignore[assignment]

    class DataLoader:  # type: ignore[override]
        def __init__(self, *args, **kwargs) -> None:
            raise ModuleNotFoundError(
                "ReactionDataLoader requires DeepChem/TensorFlow to be installed."
            ) from _DEEPCHEM_IMPORT_ERROR

    class SingletaskStratifiedSplitter:  # type: ignore[override]
        def __init__(self, *args, **kwargs) -> None:
            raise ModuleNotFoundError(
                "stratified_split requires DeepChem/TensorFlow to be installed."
            ) from _DEEPCHEM_IMPORT_ERROR

from deepretro.featurizers import ReactionStepFeaturizer


class ReactionDataLoader(DataLoader):
    """Load a two-column reaction CSV into a DeepChem ``DiskDataset``.

    Inherits from :class:`~deepchem.data.DataLoader` and overrides
    :meth:`_get_shards` and :meth:`_featurize_shard` to handle paired
    SMILES columns (product + reactants).  The parent's
    ``create_dataset`` method is overridden to add automatic NaN-row
    dropping (invalid SMILES) with a summary warning.

    Parameters
    ----------
    featurizer : ReactionStepFeaturizer or None, optional
        Pre-configured featurizer.  A default one (radius=2, size=2048,
        domain features on) is created when ``None``.
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
        Log a progress message every *n* shards.  Default ``1000``.

    Examples
    --------
    >>> from deepretro.data import ReactionDataLoader  # doctest: +SKIP
    >>> loader = ReactionDataLoader()  # doctest: +SKIP
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
            (all-NaN rows are excluded).
        """
        products = shard[self.product_col].tolist()
        reactants = shard[self.reactants_col].tolist()
        pairs = list(zip(products, reactants))
        X = self.featurizer.featurize(pairs)

        # Identify rows where featurization failed (all-NaN)
        nan_mask = np.isnan(X).all(axis=1)
        valid_inds = np.where(~nan_mask)[0]
        X = X[valid_inds]
        return X, valid_inds

    # create_dataset : customised to handle missing id_field & NaN warn

    def create_dataset(
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
            directory is used (cleaned up when the ``DiskDataset``
            object is garbage-collected).
        shard_size : int or None, optional
            Number of rows per shard.  Default ``1000``.

        Returns
        -------
        dataset : DiskDataset
            Disk-backed DeepChem dataset ready for splitting / training.

        Examples
        --------
        >>> loader = ReactionDataLoader()  # doctest: +SKIP
        >>> ds = loader.create_dataset("data/dataset.csv")  # doctest: +SKIP
        """
        if not isinstance(inputs, list):
            inputs = [inputs]

        total_dropped = 0

        # Generator that streams featurized shards (X, y, w, ids) into ``DiskDataset.create_dataset`` while tracking how many rows are dropped due to all-NaN features from failed featurization.
        def shard_generator():
            nonlocal total_dropped
            for shard in self._get_shards(inputs, shard_size):
                X, valid_inds = self._featurize_shard(shard)

                dropped = len(shard) - len(valid_inds)
                if dropped > 0:
                    total_dropped += dropped

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
                f"Dropped {total_dropped} rows with NaN features (invalid SMILES)."
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
    >>> import numpy as np  # doctest: +SKIP
    >>> from deepchem.data import NumpyDataset  # doctest: +SKIP
    >>> from deepretro.data import stratified_split  # doctest: +SKIP
    >>> ds = NumpyDataset(X=np.random.rand(100, 10), y=np.array([0]*50 + [1]*50).reshape(-1,1))  # doctest: +SKIP
    >>> train, valid, test = stratified_split(ds)  # doctest: +SKIP
    >>> len(train) + len(valid) + len(test) == 100  # doctest: +SKIP
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
