__all__ = [
    "extract_domain_features_single",
    "NUM_DOMAIN_FEATURES",
    "find_optimal_threshold",
]


def __getattr__(name: str):
    """Lazily import utility helpers to avoid import-time heavy dependencies."""
    if name in {"extract_domain_features_single", "NUM_DOMAIN_FEATURES"}:
        from .domain_features import (
            NUM_DOMAIN_FEATURES,
            extract_domain_features_single,
        )

        exports = {
            "extract_domain_features_single": extract_domain_features_single,
            "NUM_DOMAIN_FEATURES": NUM_DOMAIN_FEATURES,
        }
        return exports[name]
    if name == "find_optimal_threshold":
        from .metrics import find_optimal_threshold

        return find_optimal_threshold
    raise AttributeError(f"module 'deepretro.utils' has no attribute {name!r}")
