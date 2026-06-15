"""Subtype tagging utilities for hallucination-dataset annotations.

This module converts free-text ``reason`` annotations into a small,
normalized subtype taxonomy that can be used for downstream analysis
without changing the current model-training behavior.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

SUBTYPE_VALID_CLEAN = "valid_clean"
SUBTYPE_PROTECTING_GROUP_ISSUE = "protecting_group_issue"
SUBTYPE_NO_REAL_DISCONNECTION = "no_real_disconnection"
SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION = "wrong_starting_material_or_direction"
SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR = "atom_count_or_formula_error"
SUBTYPE_RING_TOPOLOGY_ERROR = "ring_topology_error"
SUBTYPE_REGIOCHEMICAL_ERROR = "regiochemical_error"
SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR = "bond_or_fragment_connectivity_error"
SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR = "unstable_intermediate_or_precursor"
SUBTYPE_UNCLASSIFIED_HALLUCINATION = "unclassified_hallucination"

ALL_SUBTYPES = (
    SUBTYPE_VALID_CLEAN,
    SUBTYPE_PROTECTING_GROUP_ISSUE,
    SUBTYPE_NO_REAL_DISCONNECTION,
    SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION,
    SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
    SUBTYPE_RING_TOPOLOGY_ERROR,
    SUBTYPE_REGIOCHEMICAL_ERROR,
    SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
    SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR,
    SUBTYPE_UNCLASSIFIED_HALLUCINATION,
)


@dataclass(frozen=True)
class SubtypeRule:
    """Map one or more reason-text markers to a normalized subtype."""

    subtype: str
    markers: tuple[str, ...]


_PRIMARY_RULES: tuple[SubtypeRule, ...] = (
    SubtypeRule(
        SUBTYPE_REGIOCHEMICAL_ERROR,
        ("position hallucination", "wrong position", "regio"),
    ),
    SubtypeRule(
        SUBTYPE_PROTECTING_GROUP_ISSUE,
        (
            "unnecessary protective group",
            "unfavourable protective group",
            "protective group",
            "protecting group",
            "methoxy",
        ),
    ),
    SubtypeRule(
        SUBTYPE_NO_REAL_DISCONNECTION,
        (
            "didn't break into smaller",
            "did not break into smaller",
            "no real progress",
            "no chemical reaction seen",
            "core unchanged",
            "pure hallucination",
        ),
    ),
    SubtypeRule(
        SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION,
        ("wrong starting material",),
    ),
    SubtypeRule(
        SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
        (
            "atom count",
            "carbon missing",
            "carbon count increased",
            "sudden change in atom count",
        ),
    ),
    SubtypeRule(
        SUBTYPE_RING_TOPOLOGY_ERROR,
        (
            "membered ring",
            "ring opening",
            "ring formation",
            "ring formed",
            "ring missing",
            "ring reduced",
        ),
    ),
    SubtypeRule(
        SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
        (
            "attachments hallucinated",
            "bond between",
            "covalency hallucination",
            "connecting carbons",
            "splits into two fragments",
            "fragment breaks into two parts",
        ),
    ),
    SubtypeRule(
        SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR,
        ("unstable",),
    ),
    SubtypeRule(
        SUBTYPE_VALID_CLEAN,
        (
            "default-clean",
            "overall chemistry looks good",
            "proper chemistry",
            "no mistakes",
        ),
    ),
)


_SECONDARY_RULES: tuple[SubtypeRule, ...] = (
    SubtypeRule(SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR, ("unstable",)),
    SubtypeRule(
        SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
        ("atom count", "carbon missing", "carbon count increased"),
    ),
    SubtypeRule(
        SUBTYPE_RING_TOPOLOGY_ERROR,
        (
            "membered ring",
            "ring opening",
            "ring formation",
            "ring missing",
            "ring reduced",
        ),
    ),
    SubtypeRule(
        SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
        ("attachments hallucinated", "bond between", "connecting carbons"),
    ),
)


def _normalize_text(reason: str | None) -> str:
    """Return a case-normalized reason string."""

    return (reason or "").strip().lower()


def _match_subtype(text: str, rules: Iterable[SubtypeRule]) -> str | None:
    """Return the first matching subtype from ``rules``."""

    for rule in rules:
        if any(marker in text for marker in rule.markers):
            return rule.subtype
    return None


def _match_secondary_subtype(
    text: str,
    rules: Iterable[SubtypeRule],
    primary: str,
) -> str | None:
    """Return the first matching subtype that differs from ``primary``."""

    for rule in rules:
        if rule.subtype == primary:
            continue
        if any(marker in text for marker in rule.markers):
            return rule.subtype
    return None


def assign_subtypes(
    label: str | int,
    reason: str | None,
) -> tuple[str, str | None]:
    """Assign normalized subtypes from a row label and free-text reason.

    Parameters
    ----------
    label : str or int
        Binary hallucination label. ``0`` is treated as clean/valid and
        ``1`` as hallucinated.
    reason : str or None
        Free-text annotation describing why the row was labeled.

    Returns
    -------
    tuple[str, str | None]
        ``(subtype_primary, subtype_secondary)``.
    """

    text = _normalize_text(reason)
    label_text = str(label).strip()

    if text:
        primary = _match_subtype(text, _PRIMARY_RULES)
        if primary is None:
            primary = (
                SUBTYPE_UNCLASSIFIED_HALLUCINATION
                if label_text == "1"
                else SUBTYPE_VALID_CLEAN
            )

        secondary = _match_secondary_subtype(text, _SECONDARY_RULES, primary)
        return primary, secondary

    if label_text == "1":
        return SUBTYPE_UNCLASSIFIED_HALLUCINATION, None
    return SUBTYPE_VALID_CLEAN, None


def enrich_annotation_row(row: dict[str, str]) -> dict[str, str]:
    """Return a copy of ``row`` with subtype columns added."""

    enriched = dict(row)
    primary, secondary = assign_subtypes(
        label=row.get("label", ""),
        reason=row.get("reason", ""),
    )
    enriched["subtype_primary"] = primary
    enriched["subtype_secondary"] = secondary or ""
    return enriched
