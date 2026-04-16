"""Model wrappers for retrosynthesis ML tasks."""

__all__ = ["HallucinationClassifier", "predict_single_reaction"]


def __getattr__(name: str):
    if name in __all__:
        from deepretro.models.hallucination_classifier import (
            HallucinationClassifier,
            predict_single_reaction,
        )
        return {"HallucinationClassifier": HallucinationClassifier,
                "predict_single_reaction": predict_single_reaction}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
