"""Unit tests for ReactionStepFeaturizer."""

import numpy as np
import pytest

from deepretro.featurizers.reactionstep import ReactionStepFeaturizer
from deepretro.utils import NUM_DOMAIN_FEATURES

ETHANOL = "CCO"
ETHANE_WATER = "CC.O"

# Real reaction from the DeepRetro dataset:
# pyrazole bromide + cyclohexanone → α-hydroxycyclohexyl pyrazole
PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"

INVALID_SMILES = "not_a_smiles!!!"


# feature_dim property


def test_feature_dim() -> None:
    feat = ReactionStepFeaturizer()
    assert feat.feature_dim == 2 * 2048

    feat_no_domain = ReactionStepFeaturizer(use_domain_features=False)
    assert feat_no_domain.feature_dim == 2 * 2048

    feat_custom = ReactionStepFeaturizer(size=1024, use_domain_features=True)
    assert feat_custom.feature_dim == 2 * 1024 + NUM_DOMAIN_FEATURES


# Single featurization (_featurize)


def test_featurize_single_valid() -> None:
    for use_domain in (True, False):
        feat = ReactionStepFeaturizer(use_domain_features=use_domain)
        result = feat._featurize((ETHANOL, ETHANE_WATER))
        assert isinstance(result, np.ndarray)
        assert result.shape == (feat.feature_dim,)
        assert result.sum() != 0.0


def test_featurize_single_invalid_returns_nan() -> None:
    feat = ReactionStepFeaturizer()
    for pair in [(INVALID_SMILES, ETHANE_WATER), (ETHANOL, INVALID_SMILES)]:
        result = feat._featurize(pair)
        assert result.shape == (feat.feature_dim,)
        assert np.all(np.isnan(result))


def test_featurize_single_distinct_reactions() -> None:
    feat = ReactionStepFeaturizer()
    v1 = feat._featurize((ETHANOL, ETHANE_WATER))
    v2 = feat._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    assert not np.array_equal(v1, v2)


# Batch featurization (featurize)


def test_featurize_batch_shape() -> None:
    feat = ReactionStepFeaturizer()
    X = feat.featurize([(ETHANOL, ETHANE_WATER)])
    assert X.shape == (1, feat.feature_dim)

    reactions = [(ETHANOL, ETHANE_WATER), (PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE)]
    X = feat.featurize(reactions)
    assert X.shape == (2, feat.feature_dim)


def test_featurize_batch_matches_single() -> None:
    feat = ReactionStepFeaturizer()
    reactions = [(ETHANOL, ETHANE_WATER), (PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE)]
    X_batch = feat.featurize(reactions)
    for i, rxn in enumerate(reactions):
        np.testing.assert_array_equal(X_batch[i], feat._featurize(rxn))


def test_featurize_batch_mixed_valid_invalid() -> None:
    feat = ReactionStepFeaturizer()
    reactions = [(ETHANOL, ETHANE_WATER), (INVALID_SMILES, ETHANE_WATER)]
    X = feat.featurize(reactions)
    assert X.shape == (2, feat.feature_dim)
    assert not np.any(np.isnan(X[0]))  # valid row has no NaNs
    assert np.all(np.isnan(X[1]))  # invalid row is all NaN


# Reproducibility


def test_reproducibility() -> None:
    feat = ReactionStepFeaturizer()
    v1 = feat._featurize((ETHANOL, ETHANE_WATER))
    v2 = feat._featurize((ETHANOL, ETHANE_WATER))
    np.testing.assert_array_equal(v1, v2)

    f2 = ReactionStepFeaturizer(radius=2, use_domain_features=False)
    f3 = ReactionStepFeaturizer(radius=3, use_domain_features=False)
    v2 = f2._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    v3 = f3._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    assert not np.array_equal(v2, v3)


def test_nested_base_featurizer_rejected() -> None:
    import pytest

    with pytest.raises(ValueError, match="cannot be used as a base"):
        ReactionStepFeaturizer(base_featurizer="reactionfeaturizer")


def test_sparse_fingerprints_rejected() -> None:
    import pytest

    with pytest.raises(ValueError, match="dense fingerprints"):
        ReactionStepFeaturizer(sparse=True)


def test_graph_features_and_dimension() -> None:
    import pytest
    from deepchem.feat.graph_data import GraphData

    featurizer = ReactionStepFeaturizer(base_featurizer="molgraphconv")
    result = featurizer.featurize([("CCO", "CC.O")])
    assert isinstance(result[0], GraphData)
    assert result[0].num_nodes == 6
    with pytest.raises(ValueError, match="fixed vector dimension"):
        _ = featurizer.feature_dim


def test_graph_batch_preserves_failure_positions() -> None:
    from deepchem.feat.graph_data import GraphData

    feat = ReactionStepFeaturizer(base_featurizer="molgraphconv")
    result = feat.featurize([("CCO", "CC.O"), (INVALID_SMILES, "CC.O")])
    assert result.shape == (2,)
    assert isinstance(result[0], GraphData)
    assert result[1].size == 0


def test_empty_batch_retains_vector_width() -> None:
    feat = ReactionStepFeaturizer(size=16)
    assert feat.featurize([]).shape == (0, 32)


@pytest.mark.parametrize("use_domain_features", [False, True])
def test_registered_vector_featurizer_without_size(
    monkeypatch: pytest.MonkeyPatch, use_domain_features: bool
) -> None:
    """Custom DeepChem descriptor vectors retain their inferred width."""
    from deepchem.feat import RDKitDescriptors
    from deepretro.featurizers.reactionstep import FEATURIZER_MAP

    monkeypatch.setitem(
        FEATURIZER_MAP,
        "descriptors",
        {"class": RDKitDescriptors, "type": "1d_vector", "default_params": {}},
    )
    feat = ReactionStepFeaturizer(
        base_featurizer="descriptors", use_domain_features=use_domain_features
    )
    expected = np.concatenate(
        [feat.base_featurizer.featurize([smi])[0] for smi in (ETHANOL, ETHANE_WATER)]
    )
    result = feat.featurize([(INVALID_SMILES, ETHANE_WATER), (ETHANOL, ETHANE_WATER)])
    assert result.shape == (2, expected.size)
    assert feat.feature_dim == expected.size
    assert np.isnan(result[0]).all()
    np.testing.assert_array_equal(result[1], expected)
    assert feat.featurize([]).shape == (0, expected.size)
