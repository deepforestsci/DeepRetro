"""Unit tests for deepretro.utils.visualize."""

from __future__ import annotations

import pytest

from deepretro.utils.visualize import build_tree, mol_metadata, visualize_pathway

BENZENE = "c1ccccc1"
ETHANOL = "CCO"
ACETALDEHYDE = "CC=O"


def _result(steps: list[dict], deps: dict[str, list[str]]) -> dict:
    return {"steps": steps, "dependencies": deps}


def test_build_tree_empty_returns_none() -> None:
    assert build_tree({"steps": [], "dependencies": {}}) is None


def test_build_tree_returns_none_when_step_one_missing() -> None:
    result = {"steps": [{"step": "2", "reactants": [], "products": []}], "dependencies": {}}
    assert build_tree(result) is None


def test_build_tree_single_step() -> None:
    result = _result(
        steps=[
            {
                "step": "1",
                "reactants": [{"smiles": ETHANOL}],
                "products": [{"smiles": BENZENE}],
            },
        ],
        deps={"1": []},
    )
    root = build_tree(result)
    assert root is not None
    assert root.is_target is True
    assert root.molecules == [{"smiles": BENZENE}]
    assert len(root.children) == 1
    assert root.children[0].molecules == [{"smiles": ETHANOL}]


def test_build_tree_skips_dependency_ids_missing_from_step_map() -> None:
    result = _result(
        steps=[
            {
                "step": "1",
                "reactants": [{"smiles": ETHANOL}],
                "products": [{"smiles": ACETALDEHYDE}],
            }
        ],
        deps={"1": ["99"]},
    )
    root = build_tree(result)
    assert root is not None
    assert len(root.children) == 1
    assert root.children[0].children == []


def test_mol_metadata_uses_provided_values() -> None:
    formula, mass_str = mol_metadata(
        {
            "product_metadata": {"chemical_formula": "C2H6O", "mass": 46.07},
        }
    )
    assert formula == "C2H6O"
    assert mass_str == "46.1 g/mol"


def test_mol_metadata_falls_back_to_rdkit() -> None:
    """RDKit is used to recompute formula and mass when metadata is absent."""
    formula, mass_str = mol_metadata({"smiles": ETHANOL})
    assert formula == "C2H6O"
    assert mass_str.endswith(" g/mol")
    mass_value = float(mass_str.split()[0])
    assert mass_value == pytest.approx(46.04, abs=0.05)


def test_visualize_pathway_empty_returns_placeholder() -> None:
    img = visualize_pathway({"steps": [], "dependencies": {}})
    assert img.size == (400, 200)
    assert img.mode == "RGB"


def test_visualize_pathway_smoke() -> None:
    result = _result(
        steps=[
            {
                "step": "1",
                "reactants": [{"smiles": ETHANOL}],
                "products": [{"smiles": BENZENE}],
            },
        ],
        deps={"1": []},
    )
    img = visualize_pathway(result)
    assert img.mode == "RGB"
    assert img.size[0] > 0 and img.size[1] > 0


def test_visualize_pathway_multistep_smoke() -> None:
    result = _result(
        steps=[
            {
                "step": "1",
                "reactants": [{"smiles": ETHANOL}],
                "products": [{"smiles": BENZENE}],
            },
            {
                "step": "2",
                "reactants": [{"smiles": "CO"}],
                "products": [{"smiles": ETHANOL}],
            },
        ],
        deps={"1": ["2"], "2": []},
    )
    img = visualize_pathway(result)
    assert img.mode == "RGB"
    assert img.size[0] >= 400
    assert img.size[1] > 0
