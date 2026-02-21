"""Unit tests for ReactionStepFeaturizer."""

import numpy as np
import pytest

from deepretro.featurizers.reactionstep import ReactionStepFeaturizer
from deepretro.utils.feat_extract import NUM_DOMAIN_FEATURES

# Fixtures

ETHANOL = "CCO"
ETHANE_WATER = "CC.O"

BENZENE = "c1ccccc1"
CHLOROBENZENE = "c1ccccc1Cl"

INVALID_SMILES = "not_a_smiles!!!"


@pytest.fixture
def feat():
    """Default featurizer (radius=2, size=2048, domain features on)."""
    return ReactionStepFeaturizer()


@pytest.fixture
def feat_no_domain():
    """Featurizer without domain features."""
    return ReactionStepFeaturizer(use_domain_features=False)


# feature_dim property


class TestFeatureDim:
    def test_with_domain_features(self, feat):
        assert feat.feature_dim == 2 * 2048 + NUM_DOMAIN_FEATURES

    def test_without_domain_features(self, feat_no_domain):
        assert feat_no_domain.feature_dim == 2 * 2048

    def test_custom_size(self):
        f = ReactionStepFeaturizer(size=1024, use_domain_features=True)
        assert f.feature_dim == 2 * 1024 + NUM_DOMAIN_FEATURES


# Single featurization (_featurize)


class TestFeaturizeSingle:
    def test_returns_numpy_array(self, feat):
        result = feat._featurize((ETHANOL, ETHANE_WATER))
        assert isinstance(result, np.ndarray)

    def test_correct_shape(self, feat):
        result = feat._featurize((ETHANOL, ETHANE_WATER))
        assert result.shape == (feat.feature_dim,)

    def test_correct_shape_no_domain(self, feat_no_domain):
        result = feat_no_domain._featurize((ETHANOL, ETHANE_WATER))
        assert result.shape == (feat_no_domain.feature_dim,)

    def test_not_all_zeros_for_valid_smiles(self, feat):
        result = feat._featurize((ETHANOL, ETHANE_WATER))
        assert result.sum() != 0.0

    def test_invalid_product_returns_zeros(self, feat):
        result = feat._featurize((INVALID_SMILES, ETHANE_WATER))
        assert np.all(result == 0.0)

    def test_invalid_reactant_returns_zeros(self, feat):
        result = feat._featurize((ETHANOL, INVALID_SMILES))
        assert np.all(result == 0.0)

    def test_different_reactions_produce_different_vectors(self, feat):
        v1 = feat._featurize((ETHANOL, ETHANE_WATER))
        v2 = feat._featurize((BENZENE, CHLOROBENZENE))
        assert not np.array_equal(v1, v2)


# Batch featurization (featurize)


class TestFeaturizeBatch:
    def test_batch_shape(self, feat):
        reactions = [(ETHANOL, ETHANE_WATER), (BENZENE, CHLOROBENZENE)]
        X = feat.featurize(reactions)
        assert X.shape == (2, feat.feature_dim)

    def test_single_item_batch(self, feat):
        X = feat.featurize([(ETHANOL, ETHANE_WATER)])
        assert X.shape == (1, feat.feature_dim)

    def test_batch_matches_single(self, feat):
        reactions = [(ETHANOL, ETHANE_WATER), (BENZENE, CHLOROBENZENE)]
        X_batch = feat.featurize(reactions)
        for i, rxn in enumerate(reactions):
            np.testing.assert_array_equal(X_batch[i], feat._featurize(rxn))

    def test_mixed_valid_invalid(self, feat):
        reactions = [(ETHANOL, ETHANE_WATER), (INVALID_SMILES, ETHANE_WATER)]
        X = feat.featurize(reactions)
        assert X.shape == (2, feat.feature_dim)
        assert X[0].sum() != 0.0  # valid row is non-zero
        assert np.all(X[1] == 0.0)  # invalid row is zeros


# Reproducibility


class TestReproducibility:
    def test_same_input_same_output(self, feat):
        v1 = feat._featurize((ETHANOL, ETHANE_WATER))
        v2 = feat._featurize((ETHANOL, ETHANE_WATER))
        np.testing.assert_array_equal(v1, v2)

    def test_different_radius_different_output(self):
        f2 = ReactionStepFeaturizer(radius=2, use_domain_features=False)
        f3 = ReactionStepFeaturizer(radius=3, use_domain_features=False)
        v2 = f2._featurize((BENZENE, CHLOROBENZENE))
        v3 = f3._featurize((BENZENE, CHLOROBENZENE))
        assert not np.array_equal(v2, v3)
