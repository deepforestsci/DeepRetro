"""Model wrappers for retrosynthesis ML tasks."""

from __future__ import annotations

from deepretro.models.hallucination_helpers import (
    build_ml_checker,
    filter_with_checker,
    resolve_hallucination,
)

__all__ = [
    "HallucinationClassifier",
    "build_ml_checker",
    "filter_with_checker",
    "predict_single_reaction",
    "resolve_hallucination",
]


def __getattr__(name: str):
    # HallucinationClassifier and predict_single_reaction live in
    # hallucination_classifier which pulls in deepchem — keep that lazy.
    if name in {"HallucinationClassifier", "predict_single_reaction"}:
        from deepretro.models import hallucination_classifier

        return getattr(hallucination_classifier, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
