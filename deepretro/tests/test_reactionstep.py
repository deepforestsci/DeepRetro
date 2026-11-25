"""Unit tests for ReactionStepFeaturizer."""

import logging

import numpy as np
import pytest
from deepchem.feat import CircularFingerprint, RDKitDescriptors
from deepchem.feat.graph_data import GraphData

from deepretro.featurizers.reactionstep import FEATURIZER_MAP, ReactionStepFeaturizer

ETHANOL = "CCO"
ETHANE_WATER = "CC.O"

# Real reaction from the DeepRetro dataset:
# pyrazole bromide + cyclohexanone → α-hydroxycyclohexyl pyrazole
PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"

INVALID_SMILES = "not_a_smiles!!!"


@pytest.fixture
def descriptors_in_map(monkeypatch: pytest.MonkeyPatch) -> None:
    """Register RDKitDescriptors, a vector featurizer with no size option."""
    monkeypatch.setitem(
        FEATURIZER_MAP,
        "descriptors",
        {"class": RDKitDescriptors, "type": "1d_vector", "default_params": {}},
    )


@pytest.mark.parametrize(
    ("kwargs", "width"),
    [
        ({}, 4096),  # two 2048-bit fingerprints
        ({"size": 1024, "use_domain_features": True}, 2063),  # 2 * 1024 + 15
    ],
    ids=["default_2048_bits", "1024_bits_plus_domain"],
)
def test_feature_dim_counts_both_fingerprints(kwargs: dict, width: int) -> None:
    assert ReactionStepFeaturizer(**kwargs).feature_dim == width


def test_row_is_product_fingerprint_then_reactant_fingerprint() -> None:
    row = ReactionStepFeaturizer().featurize([(ETHANOL, ETHANE_WATER)])[0]
    ecfp = CircularFingerprint(radius=2, size=2048)
    expected = np.concatenate(ecfp.featurize([ETHANOL, ETHANE_WATER]))
    np.testing.assert_array_equal(row, expected)


def test_domain_features_end_the_row_for_ethanol_from_ethane_and_water() -> None:
    row = ReactionStepFeaturizer(use_domain_features=True).featurize(
        [(ETHANOL, ETHANE_WATER)]
    )[0]
    # C, N, O, Cl, Br deltas; bond delta -1 (CCO has 2 bonds, CC + O has 1);
    # ring, aromatic deltas; MW delta; 2 reactants; heavy-atom delta;
    # product MW, reactant MW; product rings, reactant rings.
    expected = [0, 0, 0, 0, 0, -1, 0, 0, 2.016, 2, 0, 46.069, 48.085, 0, 0]
    np.testing.assert_allclose(row[-15:], expected, atol=1e-3)


@pytest.mark.parametrize(
    "pair",
    [(INVALID_SMILES, ETHANE_WATER), (ETHANOL, INVALID_SMILES)],
    ids=["bad_product", "bad_reactants"],
)
def test_invalid_smiles_gives_all_nan_row(pair: tuple[str, str]) -> None:
    X = ReactionStepFeaturizer(size=16).featurize([pair])
    assert X.shape == (1, 32)  # 16 bits for the product + 16 for the reactants
    assert np.isnan(X).all()


def test_batch_keeps_valid_row_next_to_failed_row() -> None:
    feat = ReactionStepFeaturizer(size=16)
    X = feat.featurize([(ETHANOL, ETHANE_WATER), (INVALID_SMILES, ETHANE_WATER)])
    assert X.shape == (2, 32)  # 16 bits for the product + 16 for the reactants
    np.testing.assert_array_equal(X[0], feat.featurize([(ETHANOL, ETHANE_WATER)])[0])
    assert np.isnan(X[1]).all()


def test_batch_row_matches_featurizing_that_reaction_alone() -> None:
    feat = ReactionStepFeaturizer()
    reactions = [(ETHANOL, ETHANE_WATER), (PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE)]
    np.testing.assert_array_equal(
        feat.featurize(reactions)[1], feat.featurize(reactions[1:])[0]
    )


def test_missing_reactant_cell_gives_nan_row_and_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    # pandas reads an empty CSV cell as float NaN; DeepChem raises a TypeError on it.
    with caplog.at_level(logging.WARNING):
        X = ReactionStepFeaturizer(size=8).featurize(
            [(ETHANOL, ETHANE_WATER), (ETHANOL, float("nan"))]
        )
    assert X.shape == (2, 16)  # 8 bits for the product + 8 for the reactants
    assert np.isnan(X[1]).all()
    assert "Failed to featurize reaction 1" in caplog.text


def test_radius_changes_features() -> None:
    reaction = [(PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE)]
    radius_2 = ReactionStepFeaturizer(radius=2).featurize(reaction)
    radius_3 = ReactionStepFeaturizer(radius=3).featurize(reaction)
    assert not np.array_equal(radius_2, radius_3)


def test_empty_batch_keeps_vector_width() -> None:
    # 16 bits for the product + 16 for the reactants
    assert ReactionStepFeaturizer(size=16).featurize([]).shape == (0, 32)


def test_nested_base_featurizer_rejected() -> None:
    with pytest.raises(ValueError, match="cannot be used as a base"):
        ReactionStepFeaturizer(base_featurizer="reactionfeaturizer")


def test_sparse_fingerprints_rejected() -> None:
    with pytest.raises(ValueError, match="dense fingerprints"):
        ReactionStepFeaturizer(sparse=True)


def test_zero_bit_fingerprint_fails_the_probe_molecule() -> None:
    # RDKit cannot make a 0-bit fingerprint, so featurizing "C" fails.
    with pytest.raises(ValueError, match="the probe molecule failed"):
        ReactionStepFeaturizer(size=0)


def test_graph_joins_product_and_reactants_into_one_graph() -> None:
    result = ReactionStepFeaturizer(base_featurizer="molgraphconv").featurize(
        [(ETHANOL, ETHANE_WATER)]
    )
    assert result[0].num_nodes == 6  # 3 heavy atoms in CCO + 3 in CC.O


def test_graph_featurizer_has_no_feature_dim() -> None:
    feat = ReactionStepFeaturizer(base_featurizer="molgraphconv")
    with pytest.raises(ValueError, match="fixed vector dimension"):
        _ = feat.feature_dim


def test_graph_batch_keeps_failed_row_as_empty_array() -> None:
    feat = ReactionStepFeaturizer(base_featurizer="molgraphconv")
    result = feat.featurize([(ETHANOL, ETHANE_WATER), (INVALID_SMILES, ETHANE_WATER)])
    assert result.shape == (2,)
    assert isinstance(result[0], GraphData)
    assert result[1].size == 0


@pytest.mark.usefixtures("descriptors_in_map")
def test_descriptor_featurizer_width_is_twice_one_molecule() -> None:
    feat = ReactionStepFeaturizer(base_featurizer="descriptors")
    one_side = feat.base_featurizer.featurize([ETHANOL])[0].size
    assert feat.feature_dim == 2 * one_side


@pytest.mark.usefixtures("descriptors_in_map")
def test_domain_features_are_skipped_for_descriptors() -> None:
    with_domain = ReactionStepFeaturizer(
        base_featurizer="descriptors", use_domain_features=True
    )
    without_domain = ReactionStepFeaturizer(base_featurizer="descriptors")
    assert with_domain.feature_dim == without_domain.feature_dim


@pytest.mark.usefixtures("descriptors_in_map")
def test_featurizer_without_size_keeps_failed_row_full_width() -> None:
    feat = ReactionStepFeaturizer(base_featurizer="descriptors")
    X = feat.featurize([(INVALID_SMILES, ETHANE_WATER), (ETHANOL, ETHANE_WATER)])
    expected = np.concatenate(feat.base_featurizer.featurize([ETHANOL, ETHANE_WATER]))
    assert np.isnan(X[0]).all()
    np.testing.assert_array_equal(X[1], expected)
