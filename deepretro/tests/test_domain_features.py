"""Regression tests for the public reaction descriptor vector."""

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from deepretro.utils.domain_features import (
    NUM_DOMAIN_FEATURES,
    extract_domain_features_single,
)


@pytest.mark.parametrize(
    ("product", "reactants"),
    [
        ("invalid", "CC.O"),
        ("CCO", "invalid"),
        ("CCO", "CC.invalid"),
        ("", "CC.O"),
        ("CCO", ""),
        ("CCO", "CC."),
    ],
)
def test_invalid_component_marks_entire_row_missing(
    product: str, reactants: str
) -> None:
    """A partial mixture must not be silently encoded as a valid reaction."""
    result = extract_domain_features_single(product, reactants)
    assert result.shape == (NUM_DOMAIN_FEATURES,)
    assert np.isnan(result).all()


def test_feature_order_and_reactant_aggregation() -> None:
    """Pin vector semantics independently of the descriptor implementation."""
    result = extract_domain_features_single("CCO", "CC.O")
    product_mass = Descriptors.MolWt(Chem.MolFromSmiles("CCO"))
    reactant_mass = sum(Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in ["CC", "O"])
    expected = [
        0,
        0,
        0,
        0,
        0,
        -1,
        0,
        0,
        reactant_mass - product_mass,
        2,
        0,
        product_mass,
        reactant_mass,
        0,
        0,
    ]
    np.testing.assert_allclose(result, expected)


def test_reactant_order_does_not_change_features() -> None:
    """Mixture order is irrelevant to aggregate descriptors."""
    np.testing.assert_allclose(
        extract_domain_features_single("CCO", "CC.O"),
        extract_domain_features_single("CCO", "O.CC"),
    )
