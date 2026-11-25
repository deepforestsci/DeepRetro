"""Real-model AiZynthFinder tests for aspirin. They need no provider keys."""

from __future__ import annotations

from types import ModuleType

import pytest
from rdkit import Chem

pytestmark = [pytest.mark.integration, pytest.mark.slow]

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def test_real_expansion_policy_proposes_parsable_precursors_for_aspirin(
    az_module_with_models: ModuleType,
) -> None:
    proposals = az_module_with_models.az_single_step(ASPIRIN, top_k=10)

    assert len(proposals) >= 1
    for proposal in proposals:
        assert proposal["precursors"]
        assert all(
            Chem.MolFromSmiles(smi) is not None for smi in proposal["precursors"]
        )


def test_real_route_search_expands_aspirin(az_module_with_models: ModuleType) -> None:
    _, routes = az_module_with_models.run_az(ASPIRIN)

    assert len(routes) >= 1
    for route in routes:
        assert Chem.CanonSmiles(route["smiles"]) == Chem.CanonSmiles(ASPIRIN)
    # A searched tree has a reaction under the root; the feedstock shortcut has none.
    assert routes[0]["children"][0]["type"] == "reaction"
