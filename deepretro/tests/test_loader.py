"""Unit tests for ReactionDataLoader and stratified_split."""

import warnings

import numpy as np
import pandas as pd
import pytest
from deepchem.data import DiskDataset
from deepchem.data.data_loader import DataLoader

from deepretro.data.loader import ReactionDataLoader, stratified_split

PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"

ETHANOL = "CCO"
ETHANE_WATER = "CC.O"

BENZENE = "c1ccccc1"
METHOXYBENZENE = "c1ccccc1OC"

INVALID_SMILES = "not_a_smiles!!!"

VALID_ROWS = [
    (ETHANOL, ETHANE_WATER, 0),
    (METHOXYBENZENE, BENZENE, 0),
    (PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE, 1),
    (ETHANOL, ETHANE_WATER, 0),
    (METHOXYBENZENE, BENZENE, 1),
    (ETHANOL, ETHANE_WATER, 0),
]


@pytest.fixture
def csv_path(tmp_path):
    """Write a small valid CSV and return its path."""
    df = pd.DataFrame(VALID_ROWS, columns=["product", "reactants", "label"])
    path = tmp_path / "reactions.csv"
    df.to_csv(path, index=False)
    return str(path)


@pytest.fixture
def csv_with_invalid(tmp_path):
    """CSV that includes one row with invalid SMILES."""
    rows = VALID_ROWS + [(INVALID_SMILES, ETHANE_WATER, 0)]
    df = pd.DataFrame(rows, columns=["product", "reactants", "label"])
    path = tmp_path / "reactions_invalid.csv"
    df.to_csv(path, index=False)
    return str(path)


# ReactionDataLoader inherits from DataLoader


def test_loader_is_dataloader_subclass():
    loader = ReactionDataLoader()
    assert isinstance(loader, DataLoader)


# ReactionDataLoader.create_dataset


def test_create_dataset_returns_disk_dataset(csv_path):
    loader = ReactionDataLoader()
    ds = loader.create_dataset(csv_path)
    assert isinstance(ds, DiskDataset)


def test_create_dataset_shapes(csv_path):
    loader = ReactionDataLoader()
    ds = loader.create_dataset(csv_path)
    assert ds.X.shape == (len(VALID_ROWS), loader.featurizer.feature_dim)
    assert ds.y.shape == (len(VALID_ROWS), 1)


def test_create_dataset_custom_data_dir(csv_path, tmp_path):
    out_dir = str(tmp_path / "shards")
    loader = ReactionDataLoader()
    ds = loader.create_dataset(csv_path, data_dir=out_dir)
    assert isinstance(ds, DiskDataset)
    assert ds.data_dir == out_dir


def test_create_dataset_drops_invalid_smiles(csv_with_invalid):
    loader = ReactionDataLoader()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ds = loader.create_dataset(csv_with_invalid)
        assert len(w) == 1
        assert "Dropped 1 rows" in str(w[0].message)
    assert len(ds) == len(VALID_ROWS)


def test_create_dataset_missing_column_raises(csv_path):
    loader = ReactionDataLoader(product_col="nonexistent")
    with pytest.raises(KeyError):
        loader.create_dataset(csv_path)


def test_create_dataset_with_custom_shard_size(csv_path):
    """shard_size passed to create_dataset controls chunking."""
    loader = ReactionDataLoader()
    ds = loader.create_dataset(csv_path, shard_size=2)
    assert isinstance(ds, DiskDataset)
    assert ds.get_number_shards() >= 2
    assert len(ds) == len(VALID_ROWS)


# stratified_split


@pytest.fixture
def large_csv_path(tmp_path):
    """Write a balanced 40-row CSV (20 per class) for split tests."""
    rows = (
        [(ETHANOL, ETHANE_WATER, 0)] * 10
        + [(METHOXYBENZENE, BENZENE, 0)] * 10
        + [(PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE, 1)] * 10
        + [(METHOXYBENZENE, BENZENE, 1)] * 10
    )
    df = pd.DataFrame(rows, columns=["product", "reactants", "label"])
    path = tmp_path / "reactions_large.csv"
    df.to_csv(path, index=False)
    return str(path)


def test_split_preserves_total_count_and_shards(large_csv_path):
    """stratified_split should not drop rows and works with DiskDataset shards."""
    loader = ReactionDataLoader()
    ds = loader.create_dataset(large_csv_path, shard_size=8)
    assert isinstance(ds, DiskDataset)
    assert ds.get_number_shards() >= 2

    train, valid, test = stratified_split(ds)
    assert isinstance(train, DiskDataset)
    assert isinstance(valid, DiskDataset)
    assert isinstance(test, DiskDataset)
    assert len(train) + len(valid) + len(test) == len(ds)
