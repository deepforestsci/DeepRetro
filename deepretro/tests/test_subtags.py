"""Tests for hallucination subtag taxonomy and bootstrap seed helpers."""

from deepretro.data.subtags import (
    bootstrap_seed_rows,
    classify_primary_and_secondary_tag,
    summarize_label_counts_by_field,
)


def test_classify_expert_reason_reagent_role_error() -> None:
    primary, secondary = classify_primary_and_secondary_tag(
        label=1,
        reason=(
            "Et3N treated as nitrogen source for primary amide. "
            "Bases cannot deliver an NH2; mechanistically impossible."
        ),
    )
    assert primary == "mechanistic_impossibility"
    assert secondary == "reagent_role_error"


def test_classify_expert_reason_non_progress() -> None:
    primary, secondary = classify_primary_and_secondary_tag(
        label=1,
        reason="Molecule core unchanged; only reduction, no real progress",
    )
    assert primary == "non_disconnection_or_no_progress"
    assert secondary == "core_unchanged"


def test_bootstrap_seed_rows_prefers_annotated_expert_positives_and_adds_clean_rows() -> None:
    rows = [
        {
            "product": "P1",
            "reactants": "R1",
            "label": "1",
            "target": "T1",
            "source": "expert",
            "reason": "Unnecessary protective group manipulation",
        },
        {
            "product": "P2",
            "reactants": "R2",
            "label": "1",
            "target": "T2",
            "source": "expert",
            "reason": "Ring size hallucination; 7-membered ring formed",
        },
        {
            "product": "P3",
            "reactants": "R3",
            "label": "0",
            "target": "T3",
            "source": "expert",
            "reason": "",
        },
        {
            "product": "P4",
            "reactants": "R4",
            "label": "0",
            "target": "T4",
            "source": "expert",
            "reason": "",
        },
    ]

    seed = bootstrap_seed_rows(rows, max_positive=2, max_clean=1)

    assert len(seed) == 3
    assert sum(int(row["label"]) for row in seed) == 2
    assert seed[0]["primary_hallucination_tag"] == "protecting_group_strategy_error"
    assert seed[0]["secondary_hallucination_tag"] == "unnecessary_protection"
    assert seed[1]["primary_hallucination_tag"] == "structural_topology_error"
    assert seed[1]["secondary_hallucination_tag"] == "ring_size_change"
    assert seed[2]["primary_hallucination_tag"] == "clean_or_plausible_step"
    assert seed[2]["secondary_hallucination_tag"] == ""


def test_summarize_label_counts_by_field_groups_rows() -> None:
    rows = [
        {
            "label": "1",
            "source": "expert",
            "primary_hallucination_tag": "mechanistic_impossibility",
        },
        {
            "label": "0",
            "source": "expert",
            "primary_hallucination_tag": "clean_or_plausible_step",
        },
        {
            "label": "1",
            "source": "patent",
            "primary_hallucination_tag": "",
        },
    ]

    summary = summarize_label_counts_by_field(rows, "source")

    assert summary == {
        "expert": {"0": 1, "1": 1},
        "patent": {"1": 1},
    }
