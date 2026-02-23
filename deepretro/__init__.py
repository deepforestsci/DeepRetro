"""
deepretro — retrosynthesis ML utilities.

Provides DeepChem-compatible featurizers, dataset loaders, algorithms,
and model wrappers for reaction-step data.
"""

__all__ = ["ReactionStepFeaturizer", "configure_logging", "HallucinationClassifier"]

def __getattr__(name: str):
    if name == "ReactionStepFeaturizer":
        from deepretro.featurizers.reactionstep import ReactionStepFeaturizer
        return ReactionStepFeaturizer
    if name == "configure_logging":
        from deepretro.logging import configure_logging
        return configure_logging
    if name == "HallucinationClassifier":
        from deepretro.models.hallucination_classifier import HallucinationClassifier
        return HallucinationClassifier
    raise AttributeError(f"module 'deepretro' has no attribute {name!r}")