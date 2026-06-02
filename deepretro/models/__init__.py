"""Model wrappers for retrosynthesis ML tasks."""

from __future__ import annotations

__all__ = [
    "HallucinationClassifier",
    "build_ml_checker",
    "predict_single_reaction",
    "resolve_hallucination",
]


def __getattr__(name: str):
    if name in {"HallucinationClassifier", "predict_single_reaction"}:
        from deepretro.models import hallucination_classifier

        return getattr(hallucination_classifier, name)
    if name in {"build_ml_checker", "resolve_hallucination"}:
        from deepretro.models import hallucination_helpers

        return getattr(hallucination_helpers, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
