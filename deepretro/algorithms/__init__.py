"""Recursive retrosynthesis algorithms.

Combines AiZynthFinder (template-based) with LLM-based disconnection proposals
for hybrid retrosynthesis planning.

Exports
-------
rec_run_DeepRetro : Full recursive decomposition of a molecule.
single_run_DeepRetro : Single-step retrosynthesis (AZ + LLM).
MolNode, ReactionNode, ReactionMetadata : Tree node types.
ResultPair : Type alias for ``(MolNode, bool)``.
"""

from deepretro.algorithms.recursive_algo import (
    MolNode,
    ReactionMetadata,
    ReactionNode,
    ResultPair,
    rec_run_DeepRetro,
    single_run_DeepRetro,
)

__all__ = [
    "MolNode",
    "ReactionMetadata",
    "ReactionNode",
    "ResultPair",
    "rec_run_DeepRetro",
    "single_run_DeepRetro",
]
