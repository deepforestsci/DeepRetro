"""Model wrappers for retrosynthesis ML tasks."""

from __future__ import annotations

from typing import Any

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


def __getattr__(name: str) -> Any:
    """Lazily resolve deepchem-backed attributes (PEP 562 module hook).

    ``HallucinationClassifier`` and ``predict_single_reaction`` live in
    :mod:`deepretro.models.hallucination_classifier`, which imports deepchem.
    Loading them lazily keeps ``import deepretro.models`` cheap for callers
    that only need the lightweight hallucination helpers.

    Parameters
    ----------
    name : str
        Attribute name requested on the module.

    Returns
    -------
    Any
        The resolved attribute from ``hallucination_classifier``.

    Raises
    ------
    AttributeError
        If *name* is not a lazily exported attribute.

    Examples
    --------
    >>> from deepretro import models
    >>> models.HallucinationClassifier  # doctest: +SKIP
    <class 'deepretro.models.hallucination_classifier.HallucinationClassifier'>
    """
    # HallucinationClassifier and predict_single_reaction live in
    # hallucination_classifier which pulls in deepchem — keep that lazy.
    if name in {"HallucinationClassifier", "predict_single_reaction"}:
        from deepretro.models import hallucination_classifier

        return getattr(hallucination_classifier, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
