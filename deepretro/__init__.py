"""
deepretro — retrosynthesis ML utilities.

Provides DeepChem-compatible featurizers for reaction-step data.
"""

__all__ = ["ReactionStepFeaturizer", "configure_logging"]


def __getattr__(name: str):
    if name == "ReactionStepFeaturizer":
        from deepretro.featurizers.reactionstep import ReactionStepFeaturizer

        return ReactionStepFeaturizer
    if name == "configure_logging":
        from deepretro.logging import configure_logging

        return configure_logging
    raise AttributeError(f"module 'deepretro' has no attribute {name!r}")
