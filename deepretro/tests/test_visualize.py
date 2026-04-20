"""Unit tests for deepretro.utils.visualize.

Covers the parsing and end-to-end path only.  Layout math and
PIL/RDKit rendering wrappers are not unit-tested — they're either
trivial geometry or thin wrappers whose failures surface visually.
"""

from __future__ import annotations

from deepretro.utils.visualize import (
    build_tree,
    mol_metadata,
    visualize_pathway,
)


BENZENE = "c1ccccc1"
ETHANOL = "CCO"


def _result(steps: list[dict], deps: dict[str, list[str]]) -> dict:
    return {"steps": steps, "dependencies": deps}


def test_build_tree_empty_returns_none() -> None:
    assert build_tree({"steps": [], "dependencies": {}}) is None


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


def test_mol_metadata_uses_provided_values() -> None:
    formula, mass_str = mol_metadata({
        "product_metadata": {"chemical_formula": "C2H6O", "mass": 46.07},
    })
    assert formula == "C2H6O"
    assert "46" in mass_str


def test_mol_metadata_falls_back_to_rdkit() -> None:
    formula, _ = mol_metadata({"smiles": ETHANOL})
    assert formula == "C2H6O"


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
