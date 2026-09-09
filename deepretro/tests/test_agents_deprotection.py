"""Structural candidate and dispatch checks for chemical deprotection mode."""

from __future__ import annotations

from typing import Any

import pytest
from rdkit import Chem

from deepretro.agents.protection import ProtectionContext
from deepretro.agents.tools import build_tool_registry


def canonical(smiles: str) -> str:
    """Return canonical isomeric SMILES for a known-valid test molecule."""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


@pytest.mark.parametrize(
    ("group", "protected", "product"),
    [
        ("OMe", "COc1ccccc1", "Oc1ccccc1"),
        ("OEt", "CCOc1ccccc1", "Oc1ccccc1"),
        ("OBn", "CCOCc1ccccc1", "CCO"),
        ("TBS", "CCO[Si](C)(C)C(C)(C)C", "CCO"),
        ("Boc", "CCNC(=O)OC(C)(C)C", "CCN"),
        ("Boc", "CN(C)C(=O)OC(C)(C)C", "CNC"),
        ("Cbz", "CCNC(=O)OCc1ccccc1", "CCN"),
        ("Cbz", "CN(C)C(=O)OCc1ccccc1", "CNC"),
    ],
)
def test_supported_groups_produce_valid_single_site_candidates(
    group: str, protected: str, product: str
) -> None:
    result = ProtectionContext().handle_deprotection(
        protected, mode="propose", groups=[group]
    )
    assert result["mode"] == "propose"
    assert result["direction"] == "forward"
    assert len(result["candidates"]) == 1
    candidate = result["candidates"][0]
    assert candidate["protected_smiles"] == canonical(protected)
    assert candidate["product_smiles"] == canonical(product)
    assert candidate["removed_group"] == group
    assert (
        candidate["reaction_smiles"] == f"{canonical(protected)}>>{canonical(product)}"
    )
    assert candidate["conditions"] is None
    assert candidate["status"] == "structural_candidate"
    original = Chem.MolFromSmiles(candidate["protected_smiles"])
    anchor = original.GetAtomWithIdx(candidate["site_atom_index"])
    assert anchor.GetAtomicNum() == (7 if group in {"Boc", "Cbz"} else 8)
    assert anchor.GetAtomMapNum() == candidate["site_atom_map"]
    mol = Chem.MolFromSmiles(candidate["product_smiles"])
    assert mol is not None
    assert all(
        a.GetAtomicNum() != 0 and a.GetNumRadicalElectrons() == 0
        for a in mol.GetAtoms()
    )


@pytest.mark.parametrize(
    ("protected", "product", "group", "anchor_map"),
    [
        (
            "[CH3:1][NH:2][C:3](=[O:4])[O:5][C:6]([CH3:7])([CH3:8])[CH3:9]",
            "[CH3:1][NH2:2]",
            "Boc",
            2,
        ),
        (
            "[CH3:1][N:2]([CH3:10])[C:3](=[O:4])[O:5][C:6]([CH3:7])([CH3:8])[CH3:9]",
            "[CH3:1][NH:2][CH3:10]",
            "Boc",
            2,
        ),
        ("[CH3:11][O:12][c:13]1ccccc1", "[OH:12][c:13]1ccccc1", "OMe", 12),
    ],
)
def test_mapped_explicit_hydrogens_are_updated_without_radicals(
    protected: str, product: str, group: str, anchor_map: int
) -> None:
    candidate = ProtectionContext().handle_deprotection(
        protected, mode="propose", groups=[group]
    )["candidates"][0]
    assert candidate["product_smiles"] == canonical(product)
    assert candidate["site_atom_map"] == anchor_map
    mol = Chem.MolFromSmiles(candidate["product_smiles"])
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in mol.GetAtoms())


def test_proposal_preserves_scaffold_stereochemistry_and_salts() -> None:
    protected = "CO[C@@H](C)C(=O)O.[Na+].[Cl-]"
    candidate = ProtectionContext().handle_deprotection(
        protected, mode="propose", groups=["OMe"]
    )["candidates"][0]
    assert candidate["product_smiles"] == canonical("O[C@@H](C)C(=O)O.[Na+].[Cl-]")


def test_each_candidate_removes_only_one_of_two_different_groups() -> None:
    protected = "COc1ccc(OCc2ccccc2)cc1"
    candidates = ProtectionContext().handle_deprotection(protected, mode="propose")[
        "candidates"
    ]
    assert len(candidates) == 2
    assert {c["removed_group"]: c["product_smiles"] for c in candidates} == {
        "OMe": canonical("Oc1ccc(OCc2ccccc2)cc1"),
        "OBn": canonical("COc1ccc(O)cc1"),
    }
    assert len({c["site_atom_index"] for c in candidates}) == 2


def test_equivalent_sites_still_have_separate_site_candidates() -> None:
    candidates = ProtectionContext().handle_deprotection(
        "COc1ccc(OC)cc1", mode="propose", groups=["OMe"]
    )["candidates"]
    assert len(candidates) == 2
    assert len({c["site_atom_index"] for c in candidates}) == 2
    assert all(c["product_smiles"] == canonical("COc1ccc(O)cc1") for c in candidates)


def test_disconnected_identical_components_keep_other_component() -> None:
    candidates = ProtectionContext().handle_deprotection(
        "COc1ccccc1.COc1ccccc1", mode="propose", groups=["OMe"]
    )["candidates"]
    assert len(candidates) == 2
    assert all(
        c["product_smiles"] == canonical("Oc1ccccc1.COc1ccccc1") for c in candidates
    )


@pytest.mark.parametrize(
    ("smiles", "groups"),
    [
        ("CCO", None),
        ("CCN", None),
        ("Oc1ccccc1", None),
        ("COc1ccccc1", []),
        ("COc1ccccc1", ["Boc"]),
    ],
)
def test_no_matching_group_returns_no_candidates(
    smiles: str, groups: list[str] | None
) -> None:
    result = ProtectionContext().handle_deprotection(
        smiles, mode="propose", groups=groups
    )
    assert result["candidates"] == []


@pytest.mark.parametrize(
    "arguments",
    [
        {"smiles": "COc1ccccc1", "mode": "invalid"},
        {"smiles": "COc1ccccc1", "mode": "propose", "groups": ["unknown"]},
        {"smiles": "COc1ccccc1", "mode": "propose", "groups": "OMe"},
        {"smiles": "[*:1]Oc1ccccc1", "mode": "propose"},
        {"smiles": "COc1ccccc1", "mode": "propose", "mask_id": "mask_any"},
        {"smiles": "invalid_smiles", "mode": "propose"},
        {"smiles": "", "mode": "propose"},
    ],
)
def test_dispatcher_rejects_invalid_proposal_requests(
    arguments: dict[str, Any],
) -> None:
    assert "error" in build_tool_registry().execute("handle_deprotection", arguments)


def test_default_restore_remains_backward_compatible_and_rejects_group_filter() -> None:
    registry = build_tool_registry()
    masked = registry.execute("handle_protection", {"smiles": "COc1ccccc1"})
    arguments = {"smiles": masked["masked_smiles"], "mask_id": masked["mask_id"]}
    restored = registry.execute("handle_deprotection", arguments)
    assert restored["smiles"] == "COc1ccccc1"
    assert restored["restored_groups"] == 1
    assert "error" in registry.execute(
        "handle_deprotection", {**arguments, "groups": ["OMe"]}
    )


def test_proposal_is_exposed_in_tool_schema_and_dispatcher() -> None:
    registry = build_tool_registry()
    schema = next(
        s["function"]
        for s in registry.schemas
        if s["function"]["name"] == "handle_deprotection"
    )
    assert "mode" in schema["parameters"]["properties"]
    assert "mask_id" not in schema["parameters"]["required"]
    result = registry.execute(
        "handle_deprotection", {"smiles": "COc1ccccc1", "mode": "propose"}
    )
    assert result["candidates"][0]["product_smiles"] == "Oc1ccccc1"


def test_proposing_does_not_change_an_existing_restoration_handle() -> None:
    registry = build_tool_registry()
    masked = registry.execute("handle_protection", {"smiles": "COc1ccccc1"})
    proposed = registry.execute(
        "handle_deprotection", {"smiles": "CCNC(=O)OC(C)(C)C", "mode": "propose"}
    )
    assert proposed["candidates"][0]["product_smiles"] == "CCN"
    restored = registry.execute(
        "handle_deprotection",
        {"smiles": masked["masked_smiles"], "mask_id": masked["mask_id"]},
    )
    assert restored["smiles"] == masked["original_smiles"]
