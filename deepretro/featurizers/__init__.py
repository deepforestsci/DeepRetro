"""Featurizers for reaction-step data.

Exports :class:`ReactionStepFeaturizer` for converting product-reactant pairs
into fixed-size vectors, and :data:`NUM_DOMAIN_FEATURES` (15) for the number
of hand-crafted domain features appended when ``use_domain_features=True``.
"""

from deepretro.featurizers.reactionstep import ReactionStepFeaturizer

__all__ = ["ReactionStepFeaturizer"]
