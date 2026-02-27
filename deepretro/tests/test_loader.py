"""Unit tests for ReactionDataLoader and split_dataset."""

import warnings

import numpy as np
import pandas as pd
import pytest
from deepchem.data import NumpyDataset

from deepretro.data.loader import ReactionDataLoader, split_dataset

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


# ReactionDataLoader.create_dataset


def test_create_dataset_shapes(csv_path):
    loader = ReactionDataLoader()
    ds = loader.create_dataset(csv_path)
    assert isinstance(ds, NumpyDataset)
    assert ds.X.shape == (len(VALID_ROWS), loader.featurizer.feature_dim)
    assert ds.y.shape == (len(VALID_ROWS), 1)


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


# split_dataset


def test_split_preserves_total_count():
    ds = NumpyDataset(
        X=np.random.RandomState(0).rand(100, 10),
        y=np.array([0] * 50 + [1] * 50).reshape(-1, 1),
    )
    train, valid, test = split_dataset(ds)
    assert isinstance(train, NumpyDataset)
    assert len(train) + len(valid) + len(test) == len(ds)
