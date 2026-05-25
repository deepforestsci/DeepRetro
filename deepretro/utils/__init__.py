from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .domain_features import extract_domain_features_single, NUM_DOMAIN_FEATURES
    from .metrics import find_optimal_threshold
else:
    try:
        from .domain_features import extract_domain_features_single, NUM_DOMAIN_FEATURES
    except Exception:
        pass

    try:
        from .metrics import find_optimal_threshold
    except Exception:
        pass

__all__ = [
    "extract_domain_features_single",
    "NUM_DOMAIN_FEATURES",
    "find_optimal_threshold",
    "visualize_pathway",
]


def __getattr__(name: str):
    if name == "visualize_pathway":
        from deepretro.utils.visualize import visualize_pathway

        return visualize_pathway
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
