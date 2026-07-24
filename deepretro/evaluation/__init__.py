"""Offline plausibility and diversity evaluation for DeepRetro predictions.

This package deliberately does not participate in the production ``src``
retrosynthesis runtime.  It supplies reproducible, research-facing evaluation
utilities for ranked single-step precursor predictions.
"""

from deepretro.evaluation.canonical import canonicalize_prediction
from deepretro.evaluation.cache import ScoreCache
from deepretro.evaluation.aggregate import aggregate_benchmark, evaluate_target
from deepretro.evaluation.io import load_predictions, load_scores, load_targets, write_scores
from deepretro.evaluation.runner import score_predictions
from deepretro.evaluation.report import write_benchmark_report, write_score_provenance
from deepretro.evaluation.chemcensor import OfficialChemCensorScorer
from deepretro.evaluation.types import (
    CanonicalPrediction,
    ChemCensorProvenance,
    ChemCensorScore,
    EvaluationTarget,
    SingleStepPrediction,
)

__all__ = [
    "CanonicalPrediction",
    "ChemCensorProvenance",
    "ChemCensorScore",
    "EvaluationTarget",
    "OfficialChemCensorScorer",
    "ScoreCache",
    "aggregate_benchmark",
    "SingleStepPrediction",
    "canonicalize_prediction",
    "evaluate_target",
    "load_predictions",
    "load_scores",
    "load_targets",
    "score_predictions",
    "write_benchmark_report",
    "write_score_provenance",
    "write_scores",
]
