"""
deepretro — retrosynthesis ML utilities.

Provides DeepChem-compatible featurizers for reaction-step data, recursive
planning algorithms, and LLM/AiZynthFinder integration.

Subpackages
-----------
algorithms : Recursive retrosynthesis (AiZynthFinder + LLM)
featurizers : Reaction-step featurizers for ML models
utils : LLM calls, AiZynthFinder, and domain feature extraction
"""

from deepretro.algorithms import (
    rec_run_DeepRetro,
    single_run_DeepRetro,
)
from deepretro.featurizers.reactionstep import ReactionStepFeaturizer
from deepretro.utils import (
    llm_pipeline,
    run_az,
)

__all__ = [
    "ReactionStepFeaturizer",
    "rec_run_DeepRetro",
    "single_run_DeepRetro",
    "llm_pipeline",
    "run_az",
]
