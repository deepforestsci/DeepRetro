"""Unit tests for deepretro.utils.visualize."""

from __future__ import annotations

import pytest

import deepretro.utils.visualize as visualize
from deepretro.utils.visualize import build_tree, mol_metadata, render_mol, visualize_pathway

BENZENE = "c1ccccc1"
ETHANOL = "CCO"
ACETALDEHYDE = "CC=O"


def _result(steps: list[dict], deps: dict[str, list[str]]) -> dict:
    """Build a minimal viewer-ready pathway payload for test cases."""
    return {"steps": steps, "dependencies": deps}


def test_build_tree_empty_returns_none() -> None:
    """Empty pathway payloads should not produce a visualization tree."""
    assert build_tree({"steps": [], "dependencies": {}}) is None


def test_build_tree_returns_none_when_step_one_missing() -> None:
    """The viewer contract requires step ``1`` as the root reaction node."""
    result = {"steps": [{"step": "2", "reactants": [], "products": []}], "dependencies": {}}
    assert build_tree(result) is None


def test_build_tree_single_step() -> None:
    """A one-step route should render the target plus one reactant node."""
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
    """Dangling dependency ids should be ignored rather than crashing tree build."""
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


def test_build_tree_uses_product_for_terminal_product_only_step() -> None:
    """Terminal product-only steps should render their product instead of disappearing."""
    result = _result(
        steps=[
            {
                "step": "1",
                "reactants": [{"smiles": ETHANOL}],
                "products": [{"smiles": BENZENE}],
            },
            {
                "step": "2",
                "reactants": [],
                "products": [{"smiles": ETHANOL}],
            },
        ],
        deps={"1": ["2"], "2": []},
    )
    root = build_tree(result)

    assert root is not None
    assert len(root.children) == 1
    assert len(root.children[0].children) == 1
    assert root.children[0].children[0].molecules == [{"smiles": ETHANOL}]


def test_build_tree_keeps_nonterminal_product_only_step_empty() -> None:
    """Only leaf product-only steps should fall back to products for rendering."""
    result = _result(
        steps=[
            {
                "step": "1",
                "reactants": [{"smiles": BENZENE}],
                "products": [{"smiles": ACETALDEHYDE}],
            },
            {
                "step": "2",
                "reactants": [],
                "products": [{"smiles": ETHANOL}],
            },
            {
                "step": "3",
                "reactants": [{"smiles": "CO"}],
                "products": [{"smiles": ETHANOL}],
            },
        ],
        deps={"1": ["2"], "2": ["3"], "3": []},
    )
    root = build_tree(result)

    assert root is not None
    assert root.children[0].children[0].molecules == []


def test_mol_metadata_uses_provided_values() -> None:
    """Embedded metadata should be preferred over recomputation."""
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


def test_mol_metadata_returns_smiles_for_invalid_smiles() -> None:
    """Invalid SMILES should degrade to a readable fallback instead of raising."""
    formula, mass_str = mol_metadata({"smiles": "not-a-smiles"})

    assert formula == "not-a-smiles"
    assert mass_str == "?"


def test_visualize_pathway_requires_pillow(monkeypatch: pytest.MonkeyPatch) -> None:
    """Visualization entrypoints should raise a clear error when Pillow is missing."""
    monkeypatch.setattr(visualize, "Image", None)
    monkeypatch.setattr(visualize, "ImageDraw", None)
    monkeypatch.setattr(visualize, "ImageFont", None)

    with pytest.raises(ImportError, match="requires Pillow"):
        visualize_pathway({"steps": [], "dependencies": {}})


def test_resample_filter_falls_back_without_resampling(monkeypatch: pytest.MonkeyPatch) -> None:
    """Older Pillow releases should still expose a usable resize filter."""

    class LegacyImageModule:
        LANCZOS = 7

    monkeypatch.setattr(visualize, "Image", LegacyImageModule)

    assert visualize.resample_filter() == 7


def test_visualize_pathway_empty_returns_placeholder() -> None:
    """Empty routes should produce a small placeholder image."""
    img = visualize_pathway({"steps": [], "dependencies": {}})
    assert img.size == (400, 200)
    assert img.mode == "RGB"


def test_visualize_pathway_missing_step_one_returns_placeholder() -> None:
    """Routes without the required root step should fall back to the placeholder image."""
    img = visualize_pathway({"steps": [{"step": "2", "reactants": [], "products": []}], "dependencies": {}})

    assert img.size == (400, 200)
    assert img.mode == "RGB"


def test_visualize_pathway_smoke() -> None:
    """A simple one-step route should render a non-empty RGB image."""
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


def test_visualize_pathway_multistep_uses_roomier_default_canvas() -> None:
    """Multi-step routes should now render with a roomier default canvas."""
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
    assert img.size[0] >= 1200
    assert img.size[1] >= 450


def test_resample_filter_prefers_modern_resampling_api(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Modern Pillow-style modules should prefer ``Image.Resampling.LANCZOS``."""

    class ModernResampling:
        LANCZOS = 9

    class ModernImageModule:
        Resampling = ModernResampling
        LANCZOS = 7

    monkeypatch.setattr(visualize, "Image", ModernImageModule)

    assert visualize.resample_filter() == 9


def test_get_font_falls_back_to_default_font(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Font loading should degrade to Pillow's default when TrueType fonts fail."""

    class FakeImageFontModule:
        @staticmethod
        def truetype(_name: str, _size: int) -> object:
            raise OSError("font missing")

        @staticmethod
        def load_default() -> object:
            return "default-font"

    monkeypatch.setattr(visualize, "Image", object())
    monkeypatch.setattr(visualize, "ImageFont", FakeImageFontModule)

    assert visualize.get_font(14) == "default-font"


def test_render_mol_returns_none_for_invalid_smiles() -> None:
    """Invalid SMILES should skip structure rendering instead of raising."""
    assert render_mol("not-a-smiles", 128) is None
