"""Unit tests for ReactionDataLoader and stratified_split."""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from deepchem.data import NumpyDataset
from deepchem.feat.graph_data import GraphData

from deepretro.data.loader import ReactionDataLoader, stratified_split
from deepretro.featurizers import ReactionStepFeaturizer

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
def csv_path(tmp_path: Path) -> str:
    """Write a small valid CSV and return its path."""
    df = pd.DataFrame(VALID_ROWS, columns=["product", "reactants", "label"])
    path = tmp_path / "reactions.csv"
    df.to_csv(path, index=False)
    return str(path)


@pytest.fixture
def csv_with_invalid(tmp_path: Path) -> str:
    """CSV with two failing rows (bad product, bad reactants) in the middle."""
    bad_rows = [(INVALID_SMILES, ETHANE_WATER, 1), (ETHANOL, INVALID_SMILES, 1)]
    rows = VALID_ROWS[:2] + bad_rows + VALID_ROWS[2:]
    df = pd.DataFrame(rows, columns=["product", "reactants", "label"])
    path = tmp_path / "reactions_invalid.csv"
    df.to_csv(path, index=False)
    return str(path)


def test_create_dataset_gives_row_ids_0_to_n_when_no_id_column(csv_path: str) -> None:
    ds = ReactionDataLoader().create_dataset(csv_path)
    np.testing.assert_array_equal(ds.ids, [0, 1, 2, 3, 4, 5])


def test_create_dataset_has_one_feature_row_per_csv_row(csv_path: str) -> None:
    ds = ReactionDataLoader().create_dataset(csv_path)
    assert ds.X.shape == (6, 4096)  # two 2048-bit fingerprints per reaction


def test_create_dataset_writes_to_given_data_dir(csv_path: str, tmp_path: Path) -> None:
    out_dir = str(tmp_path / "shards")
    ds = ReactionDataLoader().create_dataset(csv_path, data_dir=out_dir)
    assert ds.data_dir == out_dir


def test_create_dataset_drops_rows_with_invalid_smiles(csv_with_invalid: str) -> None:
    with pytest.warns(UserWarning, match="Dropped 2 rows"):
        ds = ReactionDataLoader().create_dataset(csv_with_invalid)
    np.testing.assert_array_equal(ds.y[:, 0], [0, 0, 1, 0, 1, 0])


def test_create_dataset_missing_column_raises(csv_path: str) -> None:
    loader = ReactionDataLoader(product_col="nonexistent")
    with pytest.raises(KeyError, match="nonexistent"):
        loader.create_dataset(csv_path)


def test_shard_size_splits_csv_into_shards(csv_path: str) -> None:
    ds = ReactionDataLoader().create_dataset(csv_path, shard_size=2)
    assert ds.get_number_shards() == 3  # 6 rows, 2 per shard
    assert len(ds) == 6


@pytest.fixture
def large_csv_path(tmp_path: Path) -> str:
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


def test_split_of_sharded_dataset_keeps_every_row(large_csv_path: str) -> None:
    ds = ReactionDataLoader().create_dataset(large_csv_path, shard_size=8)
    assert ds.get_number_shards() == 5  # 40 rows, 8 per shard
    train, valid, test = stratified_split(ds)
    assert len(train) + len(valid) + len(test) == 40


def test_split_puts_a_positive_in_every_split() -> None:
    # 10 positives in 100 rows; a random 80/10/10 split often leaves a split with none.
    y = np.array([1] * 10 + [0] * 90).reshape(-1, 1)
    ds = NumpyDataset(X=np.zeros((100, 1)), y=y)
    splits = stratified_split(ds, frac_train=0.8, frac_valid=0.1, frac_test=0.1)
    assert [int(s.y.sum()) for s in splits] == [8, 1, 1]


def test_loader_accepts_graph_features(csv_path: str) -> None:
    loader = ReactionDataLoader(
        featurizer=ReactionStepFeaturizer(base_featurizer="molgraphconv")
    )
    dataset = loader.create_dataset(csv_path)
    assert isinstance(dataset.X[0], GraphData)
    np.testing.assert_array_equal(dataset.y[:, 0], [0, 0, 1, 0, 1, 0])


def test_graph_loader_drops_invalid_rows(csv_with_invalid: str) -> None:
    loader = ReactionDataLoader(
        featurizer=ReactionStepFeaturizer(base_featurizer="molgraphconv")
    )
    with pytest.warns(UserWarning, match="Dropped 2 rows"):
        dataset = loader.create_dataset(csv_with_invalid)
    np.testing.assert_array_equal(dataset.y[:, 0], [0, 0, 1, 0, 1, 0])


def test_failed_first_shard_keeps_later_row_aligned(tmp_path: Path) -> None:
    path = tmp_path / "custom.csv"
    pd.DataFrame(
        {
            "product": [INVALID_SMILES, ETHANOL],
            "reactants": [ETHANE_WATER, ETHANE_WATER],
            "label": [0, 1],
            "id": ["bad", "good"],
        }
    ).to_csv(path, index=False)
    feat = ReactionStepFeaturizer(size=16)
    loader = ReactionDataLoader(featurizer=feat, id_field="id")
    with pytest.warns(UserWarning, match="Dropped 1 rows"):
        ds = loader.create_dataset(str(path), shard_size=1)
    np.testing.assert_array_equal(ds.X, feat.featurize([(ETHANOL, ETHANE_WATER)]))
    np.testing.assert_array_equal(ds.y, [[1]])
    np.testing.assert_array_equal(ds.ids, ["good"])


def test_all_rows_failed_gives_empty_dataset(tmp_path: Path) -> None:
    path = tmp_path / "custom.csv"
    pd.DataFrame(
        {"product": [INVALID_SMILES], "reactants": [ETHANE_WATER], "label": [0]}
    ).to_csv(path, index=False)
    loader = ReactionDataLoader(featurizer=ReactionStepFeaturizer(size=16))
    with pytest.warns(UserWarning, match="Dropped 1 rows"):
        ds = loader.create_dataset(str(path), shard_size=1)
    assert len(ds) == 0
    assert ds.X.size == ds.y.size == ds.w.size == ds.ids.size == 0


def test_loader_forwards_log_every_n(
    csv_path: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The loader passes its log_every_n value to the featurizer."""
    feat = ReactionStepFeaturizer(size=16)
    seen: list[int] = []
    original = feat.featurize

    def spy(pairs: Any, log_every_n: int = 1000, **kwargs: Any) -> np.ndarray:
        seen.append(log_every_n)
        return original(pairs, log_every_n=log_every_n, **kwargs)

    monkeypatch.setattr(feat, "featurize", spy)
    ReactionDataLoader(featurizer=feat, log_every_n=7).create_dataset(csv_path)
    assert seen == [7]
