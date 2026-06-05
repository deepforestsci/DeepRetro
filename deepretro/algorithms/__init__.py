"""Retrosynthesis algorithms: checkers, stability filters, and AutoSolver."""

from __future__ import annotations

from typing import Any

from .hallucination_checker import (
    calculate_hallucination_score,
    hallucination_compare_molecules,
    interpret_score,
)

from .stability_checker import (
    check_molecule_stability,
    is_valid_smiles,
    check_carbocations,
    check_carbenes,
    check_fused_cyclopentane,
)

__all__ = [
    "AutoSolver",
    "calculate_hallucination_score",
    "hallucination_compare_molecules",
    "interpret_score",
    "check_molecule_stability",
    "is_valid_smiles",
    "check_carbocations",
    "check_carbenes",
    "check_fused_cyclopentane",
]


def __getattr__(name: str) -> Any:
    """Lazily resolve :class:`AutoSolver` (PEP 562 module hook).

    ``AutoSolver`` lives in :mod:`deepretro.algorithms.autosolve`, which pulls
    in heavier retrosynthesis dependencies. Loading it lazily keeps
    ``import deepretro.algorithms`` cheap for callers that only need the
    lightweight checker helpers.

    Parameters
    ----------
    name : str
        Attribute name requested on the module.

    Returns
    -------
    Any
        The resolved attribute.

    Raises
    ------
    AttributeError
        If *name* is not a lazily exported attribute.

    Examples
    --------
    >>> from deepretro import algorithms
    >>> algorithms.AutoSolver  # doctest: +SKIP
    <class 'deepretro.algorithms.autosolve.AutoSolver'>
    """
    if name == "AutoSolver":
        from deepretro.algorithms.autosolve import AutoSolver

        return AutoSolver
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
