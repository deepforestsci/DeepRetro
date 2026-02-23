"""Unit tests for ReactionStepFeaturizer."""

import numpy as np

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

def test_feature_dim_with_domain_features():
    feat = ReactionStepFeaturizer()
    assert feat.feature_dim == 2 * 2048 + NUM_DOMAIN_FEATURES


def test_feature_dim_without_domain_features():
    feat = ReactionStepFeaturizer(use_domain_features=False)
    assert feat.feature_dim == 2 * 2048


def test_feature_dim_custom_size():
    feat = ReactionStepFeaturizer(size=1024, use_domain_features=True)
    assert feat.feature_dim == 2 * 1024 + NUM_DOMAIN_FEATURES


# Single featurization (_featurize)

def test_featurize_single_returns_numpy_array():
    feat = ReactionStepFeaturizer()
    result = feat._featurize((ETHANOL, ETHANE_WATER))
    assert isinstance(result, np.ndarray)


def test_featurize_single_correct_shape():
    feat = ReactionStepFeaturizer()
    result = feat._featurize((ETHANOL, ETHANE_WATER))
    assert result.shape == (feat.feature_dim,)


def test_featurize_single_correct_shape_no_domain():
    feat = ReactionStepFeaturizer(use_domain_features=False)
    result = feat._featurize((ETHANOL, ETHANE_WATER))
    assert result.shape == (feat.feature_dim,)


def test_featurize_single_not_all_zeros_for_valid_smiles():
    feat = ReactionStepFeaturizer()
    result = feat._featurize((ETHANOL, ETHANE_WATER))
    assert result.sum() != 0.0


def test_featurize_single_invalid_product_returns_nan():
    feat = ReactionStepFeaturizer()
    result = feat._featurize((INVALID_SMILES, ETHANE_WATER))
    assert np.all(np.isnan(result))


def test_featurize_single_invalid_reactant_returns_nan():
    feat = ReactionStepFeaturizer()
    result = feat._featurize((ETHANOL, INVALID_SMILES))
    assert np.all(np.isnan(result))


def test_featurize_single_different_reactions_produce_different_vectors():
    feat = ReactionStepFeaturizer()
    v1 = feat._featurize((ETHANOL, ETHANE_WATER))
    v2 = feat._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    assert not np.array_equal(v1, v2)


# Batch featurization (featurize)

def test_featurize_batch_shape():
    feat = ReactionStepFeaturizer()
    reactions = [(ETHANOL, ETHANE_WATER), (PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE)]
    X = feat.featurize(reactions)
    assert X.shape == (2, feat.feature_dim)


def test_featurize_batch_single_item():
    feat = ReactionStepFeaturizer()
    X = feat.featurize([(ETHANOL, ETHANE_WATER)])
    assert X.shape == (1, feat.feature_dim)


def test_featurize_batch_matches_single():
    feat = ReactionStepFeaturizer()
    reactions = [(ETHANOL, ETHANE_WATER), (PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE)]
    X_batch = feat.featurize(reactions)
    for i, rxn in enumerate(reactions):
        np.testing.assert_array_equal(X_batch[i], feat._featurize(rxn))


def test_featurize_batch_mixed_valid_invalid():
    feat = ReactionStepFeaturizer()
    reactions = [(ETHANOL, ETHANE_WATER), (INVALID_SMILES, ETHANE_WATER)]
    X = feat.featurize(reactions)
    assert X.shape == (2, feat.feature_dim)
    assert not np.any(np.isnan(X[0]))  # valid row has no NaNs
    assert np.all(np.isnan(X[1]))      # invalid row is all NaN


# Reproducibility

def test_reproducibility_same_input_same_output():
    feat = ReactionStepFeaturizer()
    v1 = feat._featurize((ETHANOL, ETHANE_WATER))
    v2 = feat._featurize((ETHANOL, ETHANE_WATER))
    np.testing.assert_array_equal(v1, v2)


def test_reproducibility_different_radius_different_output():
    f2 = ReactionStepFeaturizer(radius=2, use_domain_features=False)
    f3 = ReactionStepFeaturizer(radius=3, use_domain_features=False)
    v2 = f2._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    v3 = f3._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    assert not np.array_equal(v2, v3)
