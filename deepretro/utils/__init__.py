"""Shared reaction descriptors and metric helpers.

Exports are loaded on first access. Importing the utility namespace does not
load RDKit or NumPy. Dependency failures surface when a helper is requested.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .domain_features import NUM_DOMAIN_FEATURES, extract_domain_features_single
    from .metrics import find_optimal_threshold

__all__ = [
    "extract_domain_features_single",
    "NUM_DOMAIN_FEATURES",
    "find_optimal_threshold",
]


def __getattr__(name: str) -> Any:
    """Load a public utility when it is first requested.

    Parameters
    ----------
    name : str
        Public export name.

    Returns
    -------
    Any
        Requested function or constant.

    Raises
    ------
    AttributeError
        If the name is not a public export.
    """
    if name in {"extract_domain_features_single", "NUM_DOMAIN_FEATURES"}:
        from . import domain_features

        return getattr(domain_features, name)
    if name == "find_optimal_threshold":
        from .metrics import find_optimal_threshold

        return find_optimal_threshold
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
