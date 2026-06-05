"""Data utilities for reaction datasets and hallucination subtags."""

from __future__ import annotations

from deepretro.data.subtags import (
    PRIMARY_TO_SECONDARY,
    annotate_subtags,
    bootstrap_seed_rows,
    classify_primary_and_secondary_tag,
    summarize_label_counts_by_field,
)

__all__ = [
    "PRIMARY_TO_SECONDARY",
    "ReactionDataLoader",
    "annotate_subtags",
    "bootstrap_seed_rows",
    "classify_primary_and_secondary_tag",
    "summarize_label_counts_by_field",
    "stratified_split",
]


def __getattr__(name: str):
    if name in {"ReactionDataLoader", "stratified_split"}:
        from deepretro.data.loader import ReactionDataLoader, stratified_split

        exports = {
            "ReactionDataLoader": ReactionDataLoader,
            "stratified_split": stratified_split,
        }
        return exports[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
