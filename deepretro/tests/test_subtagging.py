"""Unit tests for hallucination subtype tagging utilities."""

import csv

from deepretro.data.subtagging import (
    SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
    SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
    SUBTYPE_NO_REAL_DISCONNECTION,
    SUBTYPE_PROTECTING_GROUP_ISSUE,
    SUBTYPE_REGIOCHEMICAL_ERROR,
    SUBTYPE_RING_TOPOLOGY_ERROR,
    SUBTYPE_UNCLASSIFIED_HALLUCINATION,
    SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR,
    SUBTYPE_VALID_CLEAN,
    SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION,
    assign_subtypes,
    enrich_annotation_row,
)
from scripts.enrich_hallucination_subtags import enrich_dataset


def test_assign_subtypes_maps_clean_reason():
    primary, secondary = assign_subtypes(
        label=0,
        reason="default-clean (post-V4 prompt successful pathway)",
    )
    assert primary == SUBTYPE_VALID_CLEAN
    assert secondary is None


def test_assign_subtypes_maps_position_hallucination():
    primary, secondary = assign_subtypes(
        label=1,
        reason="Position hallucination of functional group hydroxy",
    )
    assert primary == SUBTYPE_REGIOCHEMICAL_ERROR
    assert secondary is None


def test_assign_subtypes_maps_protecting_group_issue():
    primary, secondary = assign_subtypes(
        label=1,
        reason="Unnecessary protective group manipulation",
    )
    assert primary == SUBTYPE_PROTECTING_GROUP_ISSUE
    assert secondary is None


def test_assign_subtypes_maps_no_real_disconnection():
    primary, secondary = assign_subtypes(
        label=1,
        reason="Molecule core unchanged; only reduction, no real progress",
    )
    assert primary == SUBTYPE_NO_REAL_DISCONNECTION
    assert secondary is None


def test_assign_subtypes_maps_wrong_starting_material_with_secondary_unstable():
    primary, secondary = assign_subtypes(
        label=1,
        reason="Wrong starting material; very unstable starting material produced",
    )
    assert primary == SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION
    assert secondary == SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR


def test_assign_subtypes_maps_atom_count_and_ring_secondary():
    primary, secondary = assign_subtypes(
        label=1,
        reason="Carbon missing beside acetyl group; 6-membered N-ring reduced to 5-membered",
    )
    assert primary == SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR
    assert secondary == SUBTYPE_RING_TOPOLOGY_ERROR


def test_assign_subtypes_maps_connectivity_error():
    primary, secondary = assign_subtypes(
        label=1,
        reason="Several atom attachments hallucinated; bond between carbonyl and N",
    )
    assert primary == SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR
    assert secondary is None


def test_assign_subtypes_defaults_for_blank_reason():
    assert assign_subtypes(label=0, reason="") == (SUBTYPE_VALID_CLEAN, None)
    assert assign_subtypes(label=1, reason="") == (
        SUBTYPE_UNCLASSIFIED_HALLUCINATION,
        None,
    )


def test_enrich_annotation_row_adds_subtype_columns():
    row = {"label": "1", "reason": "Pure hallucination; molecule didn't break into smaller species"}
    enriched = enrich_annotation_row(row)
    assert enriched["subtype_primary"] == SUBTYPE_NO_REAL_DISCONNECTION
    assert enriched["subtype_secondary"] == ""


def test_enrich_dataset_writes_enriched_csv(tmp_path):
    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"
    with input_path.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames=["product", "reactants", "label", "reason"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "product": "CCO",
                "reactants": "CC.O",
                "label": "0",
                "reason": "default-clean (post-V4 prompt successful pathway)",
            }
        )
        writer.writerow(
            {
                "product": "c1ccccc1O",
                "reactants": "c1ccccc1.Br",
                "label": "1",
                "reason": "Position hallucination of functional group hydroxy",
            }
        )

    enrich_dataset(input_path, output_path)

    with output_path.open("r", newline="", encoding="utf-8") as infile:
        rows = list(csv.DictReader(infile))

    assert rows[0]["subtype_primary"] == SUBTYPE_VALID_CLEAN
    assert rows[1]["subtype_primary"] == SUBTYPE_REGIOCHEMICAL_ERROR
