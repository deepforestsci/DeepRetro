"""Regression checks for reversible protecting-group tool representations."""

from __future__ import annotations

from typing import Any

import pytest
from rdkit import Chem

from deepretro.agents.protection import ProtectionContext
from deepretro.agents.tools import build_tool_registry


def canonical(smiles: str) -> str:
    """Return canonical stereochemical SMILES for a test molecule."""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


@pytest.mark.parametrize(
    ("group", "smiles"),
    [
        ("OMe", "COc1ccccc1"),
        ("OEt", "CCOc1ccccc1"),
        ("OBn", "CCOCc1ccccc1"),
        ("Boc", "CCNC(=O)OC(C)(C)C"),
        ("Cbz", "CCNC(=O)OCc1ccccc1"),
        ("TBS", "CCO[Si](C)(C)C(C)(C)C"),
    ],
)
def test_each_motif_round_trips(group: str, smiles: str) -> None:
    context = ProtectionContext()
    masked = context.handle_protection(smiles, groups=[group])
    assert masked["matched"] is True
    assert [entry["group"] for entry in masked["groups"]] == [group]
    assert masked["masked_smiles"] != masked["original_smiles"]
    result = context.handle_deprotection(masked["masked_smiles"], masked["mask_id"])
    assert result["smiles"] == canonical(smiles)
    assert result["restored_groups"] == 1


@pytest.mark.parametrize(
    "smiles",
    [
        "COc1ccccc1.COc1ccccc1",
        "COc1ccccc1.[Na+].[Cl-]",
        "CO[C@@H](C)C(=O)O",
        "[CH3:71][O:72][C@@H:73]([CH3:74])[C:75](=[O:76])[OH:77]",
        "COc1ccc(OC)cc1",
    ],
)
def test_components_stereo_maps_and_multiplicity_survive(smiles: str) -> None:
    context = ProtectionContext()
    masked = context.handle_protection(smiles, groups=["OMe"])
    restored = context.handle_deprotection(masked["masked_smiles"], masked["mask_id"])
    assert restored["smiles"] == canonical(smiles)
    assert restored["restored_groups"] == len(masked["groups"])
    assert restored["smiles"].count(".") == smiles.count(".")


def test_default_prefers_whole_boc_over_contained_ethers() -> None:
    masked = ProtectionContext().handle_protection("CCNC(=O)OC(C)(C)C")
    assert [entry["group"] for entry in masked["groups"]] == ["Boc"]


@pytest.mark.parametrize("smiles", ["CCO", "CCN", "Oc1ccccc1", "CC(=O)O"])
def test_unprotected_functional_groups_are_not_hidden(smiles: str) -> None:
    masked = ProtectionContext().handle_protection(smiles)
    assert masked["matched"] is False
    assert masked["masked_smiles"] == canonical(smiles)


def test_explicit_empty_selection_disables_masking() -> None:
    masked = ProtectionContext().handle_protection("COc1ccccc1", groups=[])
    assert masked["matched"] is False


def test_stereogenic_hidden_group_is_left_explicit() -> None:
    """Isotope-chiral silicon must not invert when rebuilding a hidden group."""
    smiles = "C[Si@](Oc1ccccc1)([13CH3])C(C)(C)C"
    context = ProtectionContext()
    masked = context.handle_protection(smiles, groups=["TBS"])
    assert not masked["matched"]
    restored = context.handle_deprotection(masked["masked_smiles"], masked["mask_id"])
    assert restored["smiles"] == canonical(smiles)


@pytest.mark.parametrize(
    "smiles", ["", " ", "not_smiles", "CCO ethanol", "*OC", "[*:1]OC"]
)
def test_invalid_or_pre_masked_inputs_are_rejected(smiles: str) -> None:
    with pytest.raises(ValueError):
        ProtectionContext().handle_protection(smiles)


@pytest.mark.parametrize("groups", [["unknown"], "OMe", [1]])
def test_invalid_group_choices_are_rejected(groups: Any) -> None:
    with pytest.raises(ValueError):
        ProtectionContext().handle_protection("COc1ccccc1", groups=groups)


def test_unknown_and_duplicate_markers_are_rejected() -> None:
    context = ProtectionContext()
    masked = context.handle_protection("COc1ccccc1", groups=["OMe"])
    marker = masked["groups"][0]["marker"]
    for smiles in ["[*:999999]Oc1ccccc1", f"{marker}OC.CO{marker}"]:
        with pytest.raises(ValueError, match="unknown or duplicate"):
            context.handle_deprotection(smiles, masked["mask_id"])


def test_attachment_must_remain_on_original_heteroatom_type() -> None:
    context = ProtectionContext()
    masked = context.handle_protection("COc1ccccc1", groups=["OMe"])
    marker = masked["groups"][0]["marker"]
    with pytest.raises(ValueError, match="attachment"):
        context.handle_deprotection(f"CC{marker}", masked["mask_id"])


def test_restore_a_precursor_with_edited_core() -> None:
    context = ProtectionContext()
    masked = context.handle_protection("COc1ccc(C(=O)O)cc1", groups=["OMe"])
    marker = masked["groups"][0]["marker"]
    restored = context.handle_deprotection(f"{marker}Oc1ccc(CO)cc1", masked["mask_id"])
    assert restored["smiles"] == canonical("COc1ccc(CO)cc1")


def test_restore_one_precursor_subset_of_original_markers() -> None:
    context = ProtectionContext()
    masked = context.handle_protection("COc1ccc(OC)cc1", groups=["OMe"])
    assert len(masked["groups"]) == 2
    marker = masked["groups"][0]["marker"]
    restored = context.handle_deprotection(f"{marker}OCC", masked["mask_id"])
    assert restored["smiles"] == canonical("COCC")
    assert restored["restored_groups"] == 1


def test_tools_are_registered_and_round_trip_through_dispatcher() -> None:
    registry = build_tool_registry()
    masked = registry.execute(
        "handle_protection", {"smiles": "COc1ccccc1", "groups": ["OMe"]}
    )
    assert "error" not in masked
    restored = registry.execute(
        "handle_deprotection",
        {"smiles": masked["masked_smiles"], "mask_id": masked["mask_id"]},
    )
    assert restored["smiles"] == "COc1ccccc1"


def test_registry_contexts_cannot_restore_each_others_masks() -> None:
    first, second = build_tool_registry(), build_tool_registry()
    masked = first.execute(
        "handle_protection", {"smiles": "COc1ccccc1", "groups": ["OMe"]}
    )
    second.execute("handle_protection", {"smiles": "CCOCc1ccccc1", "groups": ["OBn"]})
    result = second.execute(
        "handle_deprotection",
        {"smiles": masked["masked_smiles"], "mask_id": masked["mask_id"]},
    )
    assert "error" in result
