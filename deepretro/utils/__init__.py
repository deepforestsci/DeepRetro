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
]
