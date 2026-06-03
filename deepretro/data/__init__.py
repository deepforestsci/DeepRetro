"""Data loading and splitting utilities for reaction datasets."""

__all__ = ["ReactionDataLoader", "stratified_split"]


def __getattr__(name: str):
    if name in __all__:
        from deepretro.data.loader import ReactionDataLoader, stratified_split

        exports = {
            "ReactionDataLoader": ReactionDataLoader,
            "stratified_split": stratified_split,
        }
        return exports[name]
    raise AttributeError(f"module 'deepretro.data' has no attribute {name!r}")
