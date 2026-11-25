"""Explicit real-model AiZynthFinder E2E coverage (no provider credentials)."""

from __future__ import annotations

import os
from pathlib import Path
from types import ModuleType

import pytest
from rdkit import Chem

pytestmark = [pytest.mark.integration, pytest.mark.slow]


def test_real_model_expansion_and_route_search(
    az_module: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Load real policies/stock, expand aspirin, and serialize a searched route."""
    configured = os.environ.get("AZ_TEST_MODELS_DIR")
    if not configured:
        pytest.fail("AZ E2E requires AZ_TEST_MODELS_DIR with downloaded models")
    config = Path(configured).expanduser().resolve() / "config.yml"
    assert config.is_file(), f"Missing AiZynthFinder configuration: {config}"
    assert az_module.AiZynthFinder is not None, "Install the deepretro az extra"
    monkeypatch.setattr(az_module, "AZ_MODEL_CONFIG_PATH", str(config))
    monkeypatch.setattr(az_module, "AZ_MODELS_PATH", str(config.parent))

    target = "CC(=O)Oc1ccccc1C(=O)O"
    assert not az_module.is_basic_molecule(target)
    assert target not in az_module.BASIC_MOLECULES
    proposals = az_module.az_single_step(target, top_k=10)
    assert proposals, "Real expansion policy produced no disconnections"
    for proposal in proposals:
        assert proposal["precursors"]
        assert 0 <= proposal["probability"] <= 1
        assert all(
            Chem.MolFromSmiles(smi) is not None for smi in proposal["precursors"]
        )

    solved, routes, finder = az_module._run_az_core(target)
    assert finder is not None, "E2E must not take the feedstock shortcut"
    assert isinstance(solved, bool)
    assert routes, "Search produced no serializable routes"
    for route in routes:
        assert route["type"] == "mol"
        assert Chem.CanonSmiles(route["smiles"]) == Chem.CanonSmiles(target)
