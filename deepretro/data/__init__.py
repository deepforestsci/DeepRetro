"""Data loading and annotation utilities for reaction datasets."""

__all__ = [
    "ReactionDataLoader",
    "stratified_split",
    "assign_subtypes",
    "enrich_annotation_row",
]


def __getattr__(name: str):
    """Lazily import heavy data helpers on first access."""

    if name in {"ReactionDataLoader", "stratified_split"}:
        from deepretro.data.loader import ReactionDataLoader, stratified_split

        exports = {
            "ReactionDataLoader": ReactionDataLoader,
            "stratified_split": stratified_split,
        }
        return exports[name]

    if name in {"assign_subtypes", "enrich_annotation_row"}:
        from deepretro.data.subtagging import assign_subtypes, enrich_annotation_row

        exports = {
            "assign_subtypes": assign_subtypes,
            "enrich_annotation_row": enrich_annotation_row,
        }
        return exports[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
