from .domain_features import extract_domain_features_single, NUM_DOMAIN_FEATURES
from .metrics import find_optimal_threshold

__all__ = [
    "extract_domain_features_single",
    "NUM_DOMAIN_FEATURES",
    "find_optimal_threshold",
    "autosolve",
    "autosolve_async",
]


def __getattr__(name: str):
    if name in ("autosolve", "autosolve_async"):
        from .autosolve import autosolve, autosolve_async

        globals()["autosolve"] = autosolve
        globals()["autosolve_async"] = autosolve_async
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
