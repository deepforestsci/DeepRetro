"""Unit tests for DiffFingerprint and RDKTopologicalFingerprint."""

import numpy as np

from deepretro.featurizers.rdkit_fingerprint import RDKTopologicalFingerprint
from deepretro.featurizers.diff_fingerprint import DiffFingerprint

ETHANOL = "CCO"
ETHANE_WATER = "CC.O"
PYRAZOLE_ADDUCT = "Cn1nccc1[C@]1(O)CCCC[C@H]1O"
PYRAZOLE_BROMIDE_KETONE = "Cn1nccc1Br.O=C1CCCC[C@H]1O"
INVALID = "not_a_smiles!!!"


def test_rdk_fingerprint_shape_and_determinism():
    fp = RDKTopologicalFingerprint(fpSize=1024)
    r1 = fp.featurize([ETHANOL])
    r2 = fp.featurize([ETHANOL])
    assert r1.shape == (1, 1024)
    np.testing.assert_array_equal(r1, r2)


def test_feature_dim():
    assert DiffFingerprint().feature_dim == 6 * 2048
    assert DiffFingerprint(size=1024).feature_dim == 6 * 1024


def test_valid_reaction_produces_correct_shape():
    feat = DiffFingerprint()
    result = feat._featurize((ETHANOL, ETHANE_WATER))
    assert result.shape == (feat.feature_dim,)
    assert not np.any(np.isnan(result))


def test_invalid_smiles_returns_nan():
    feat = DiffFingerprint()
    for pair in [(INVALID, ETHANE_WATER), (ETHANOL, INVALID)]:
        result = feat._featurize(pair)
        assert result.shape == (feat.feature_dim,)
        assert np.all(np.isnan(result))


def test_diff_block_equals_product_minus_reactant():
    feat = DiffFingerprint(size=512)
    block = 2 * 512

    result = feat._featurize((ETHANOL, ETHANE_WATER))

    reac = result[:block]
    prod = result[block : 2 * block]
    diff = result[2 * block : 3 * block]

    np.testing.assert_array_equal(diff, prod - reac)


def test_same_molecule_both_sides_gives_zero_diff():
    feat = DiffFingerprint(size=512)
    block = 2 * 512

    result = feat._featurize((ETHANOL, ETHANOL))
    diff = result[2 * block : 3 * block]

    np.testing.assert_array_equal(diff, np.zeros(block))


def test_diff_values_are_minus1_zero_or_plus1():
    feat = DiffFingerprint(size=2048)
    block = 2 * 2048

    result = feat._featurize((PYRAZOLE_ADDUCT, PYRAZOLE_BROMIDE_KETONE))
    diff = result[2 * block : 3 * block]

    assert np.all(np.isin(diff, [-1.0, 0.0, 1.0]))


def test_deterministic():
    feat = DiffFingerprint()
    v1 = feat._featurize((ETHANOL, ETHANE_WATER))
    v2 = feat._featurize((ETHANOL, ETHANE_WATER))
    np.testing.assert_array_equal(v1, v2)
