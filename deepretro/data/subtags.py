"""Helpers for hallucination subtag taxonomy and bootstrap seed generation."""

from __future__ import annotations

from typing import Iterable, Mapping


PRIMARY_TO_SECONDARY: dict[str, tuple[str, ...]] = {
    "structural_topology_error": (
        "ring_size_change",
        "invented_ring",
        "invented_attachment",
        "fragment_connector_hallucination",
        "major_scaffold_loss",
    ),
    "atom_or_mass_balance_error": (
        "atom_loss",
        "atom_gain",
        "hidden_deamination",
        "carbon_count_increase",
        "wrong_atom_identity",
    ),
    "mechanistic_impossibility": (
        "reagent_role_error",
        "fictional_one_pot",
        "unsupported_ring_closure",
        "forbidden_bond_rearrangement",
        "selectivity_without_basis",
    ),
    "unstable_or_nonviable_intermediate": (
        "gem_diamine_instability",
        "strained_intermediate",
        "nonisolable_precursor",
        "hydrolytically_unstable_intermediate",
    ),
    "non_disconnection_or_no_progress": (
        "core_unchanged",
        "trivial_reduction_only",
        "no_fragment_simplification",
        "overbreaking_without_strategy",
    ),
    "protecting_group_strategy_error": (
        "unnecessary_protection",
        "unnecessary_deprotection",
        "worse_protecting_group_swap",
    ),
    "regio_or_site_error": (
        "wrong_functional_group_position",
        "wrong_cleavage_site",
        "wrong_coupling_site",
    ),
    "multi_step_conflation": (
        "reduction_plus_coupling",
        "cascade_without_basis",
    ),
    "clean_or_plausible_step": (
        "reasonable_disconnection",
        "literature_like_transformation",
        "acceptable_protecting_group_step",
    ),
}


def classify_primary_and_secondary_tag(
    label: int | str,
    reason: str | None,
) -> tuple[str, str]:
    """Map a row's label/reason text onto the primary and secondary taxonomy."""
    if str(label) == "0":
        return "clean_or_plausible_step", ""

    text = (reason or "").strip().lower()
    if not text:
        return "mechanistic_impossibility", ""

    if "protective group" in text or "protecting group" in text:
        if "breaking" in text or "deprotection" in text:
            return "protecting_group_strategy_error", "unnecessary_deprotection"
        if "bigger" in text or "unfavourable protective group" in text:
            return "protecting_group_strategy_error", "worse_protecting_group_swap"
        return "protecting_group_strategy_error", "unnecessary_protection"

    if "no real progress" in text or "core unchanged" in text:
        return "non_disconnection_or_no_progress", "core_unchanged"
    if "didn't break into smaller" in text or "no chemical reaction" in text:
        return "non_disconnection_or_no_progress", "no_fragment_simplification"
    if "reduction step without significance" in text or "only reduction" in text:
        return "non_disconnection_or_no_progress", "trivial_reduction_only"
    if "unnecessary further breaking" in text:
        return "non_disconnection_or_no_progress", "overbreaking_without_strategy"

    if "et3n treated as nitrogen source" in text or "bases cannot deliver" in text:
        return "mechanistic_impossibility", "reagent_role_error"
    if "fictional one-pot" in text or "no precedented reaction" in text:
        return "mechanistic_impossibility", "fictional_one_pot"
    if "ring closure is asserted" in text or "unsupported" in text:
        return "mechanistic_impossibility", "unsupported_ring_closure"
    if "selective mono-hydrolysis" in text or "no selectivity reagent" in text:
        return "mechanistic_impossibility", "selectivity_without_basis"
    if "mechanistically impossible" in text or "cannot occur" in text:
        return "mechanistic_impossibility", ""

    if "gem-diamine" in text:
        return "unstable_or_nonviable_intermediate", "gem_diamine_instability"
    if "highly unstable intermediate" in text or "very unstable" in text:
        return "unstable_or_nonviable_intermediate", "strained_intermediate"
    if "non-isolable" in text or "not isolable" in text:
        return "unstable_or_nonviable_intermediate", "nonisolable_precursor"
    if "hydrolyses" in text or "collapses" in text:
        return "unstable_or_nonviable_intermediate", "hydrolytically_unstable_intermediate"

    if "position hallucination" in text:
        return "regio_or_site_error", "wrong_functional_group_position"
    if "should only break at" in text:
        return "regio_or_site_error", "wrong_cleavage_site"
    if "couples directly to an aryl amine site" in text:
        return "regio_or_site_error", "wrong_coupling_site"

    if "conflates two distinct transformations" in text:
        return "multi_step_conflation", "reduction_plus_coupling"
    if "at once" in text or "cascade" in text:
        return "multi_step_conflation", "cascade_without_basis"

    if (
        "ring size hallucination" in text
        or "ring formed" in text
        or "6-membered ring reduced to 5-membered" in text
    ):
        return "structural_topology_error", "ring_size_change"
    if "three-membered ring formed" in text or "ring formation" in text:
        return "structural_topology_error", "invented_ring"
    if "several atom attachments hallucinated" in text or "covalency hallucination" in text:
        return "structural_topology_error", "invented_attachment"
    if "connecting carbons" in text or "fragment breaks into two parts" in text:
        return "structural_topology_error", "fragment_connector_hallucination"
    if "whole 6-membered ring missing" in text or "molecule reduced to small starting material" in text:
        return "structural_topology_error", "major_scaffold_loss"

    if "carbon missing" in text or "nh2 disappears" in text:
        return "atom_or_mass_balance_error", "atom_loss"
    if "carbon count increased" in text or "grows by two carbons" in text:
        return "atom_or_mass_balance_error", "carbon_count_increase"
    if "methyl group hallucination" in text:
        return "atom_or_mass_balance_error", "wrong_atom_identity"
    if "silently absorbed" in text:
        return "atom_or_mass_balance_error", "hidden_deamination"

    return "mechanistic_impossibility", ""


def annotate_subtags(row: Mapping[str, str]) -> dict[str, str]:
    """Return a copy of a dataset row with taxonomy fields populated."""
    primary, secondary = classify_primary_and_secondary_tag(
        label=row.get("label", ""),
        reason=row.get("reason", ""),
    )
    enriched = dict(row)
    enriched["primary_hallucination_tag"] = primary
    enriched["secondary_hallucination_tag"] = secondary
    return enriched


def bootstrap_seed_rows(
    rows: Iterable[Mapping[str, str]],
    *,
    max_positive: int = 20,
    max_clean: int = 10,
) -> list[dict[str, str]]:
    """Build a small bootstrap seed from expert-annotated rows."""
    positive_rows: list[dict[str, str]] = []
    clean_rows: list[dict[str, str]] = []

    for row in rows:
        source = row.get("source", "")
        label = row.get("label", "")
        reason = row.get("reason", "")
        if source != "expert":
            continue
        if label == "1" and reason:
            positive_rows.append(annotate_subtags(row))
        elif label == "0":
            clean_rows.append(annotate_subtags(row))

    seed = positive_rows[:max_positive]
    seed.extend(clean_rows[:max_clean])
    return seed


def summarize_label_counts_by_field(
    rows: Iterable[Mapping[str, str]],
    field: str,
) -> dict[str, dict[str, int]]:
    """Group label counts by a metadata field."""
    summary: dict[str, dict[str, int]] = {}
    for row in rows:
        value = (row.get(field) or "").strip()
        if not value:
            continue
        label = str(row.get("label", "")).strip()
        if not label:
            continue
        label_counts = summary.setdefault(value, {})
        label_counts[label] = label_counts.get(label, 0) + 1
    return summary
