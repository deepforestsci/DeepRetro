"""Helpers for wiring hallucination checkers into the retrosynthesis pipeline.

Provides:

* :func:`build_ml_checker` — wrap a classifier into the callable
  signature the pipeline expects.
* :func:`resolve_hallucination` — turn a string to determine checker type
  into a checker callable (or ``None``) consumed by ``llm_pipeline``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any


def build_ml_checker(clf: Any) -> Callable:
    """Wrap a ``HallucinationClassifier`` into a pipeline-compatible checker.

    The returned callable has the same signature as the built-in
    heuristic checker (``(product, pathways) -> (int, list)``).
    Pathways flagged as hallucinated are dropped; if all are rejected
    the pipeline retries the LLM call.

    Parameters
    ----------
    clf : HallucinationClassifier
        A fitted classifier instance.

    Returns
    -------
    checker : Callable
        ``(product: str, pathways: list) -> tuple[int, list]``

    Examples
    --------
    >>> from deepretro.models import HallucinationClassifier  # doctest: +SKIP
    >>> clf = HallucinationClassifier()                       # doctest: +SKIP
    >>> clf.load("model_out/")                                # doctest: +SKIP
    >>> checker = build_ml_checker(clf)                       # doctest: +SKIP
    >>> status, kept = checker("CCO", [["CC", "O"]])          # doctest: +SKIP
    """
    def _checker(product: str, pathways: list) -> tuple[int, list]:
        from deepretro.utils.utils_molecule import is_valid_smiles

        valid = []
        for pathway in pathways:
            if isinstance(pathway, list):
                reactants_smi = ".".join(pathway)
            else:
                reactants_smi = pathway

            if not is_valid_smiles(reactants_smi):
                continue

            pred = clf.predict_single(product, reactants_smi)
            if not pred.get("is_hallucination", True):
                valid.append(pathway)

        return 200, valid

    return _checker


def resolve_hallucination(mode: str, classifier: str | Path | None) -> Callable | None:
    """Convert a user-facing hallucination mode into a checker callable.

    Parameters
    ----------
    mode : str
        One of ``"heuristic"``, ``"ml"``, or ``"none"``.
    classifier : str or Path or None
        Path to a saved ML model directory.  Required when *mode* is
        ``"ml"``, ignored otherwise.

    Returns
    -------
    Callable or None
        ``None`` when *mode* is ``"none"`` (skip checking).  Otherwise a
        callable with signature ``(product: str, pathways: list) -> (int, list)``
        that filters out hallucinated pathways.

    Raises
    ------
    ValueError
        If *mode* is not recognised, or *mode* is ``"ml"`` and
        *classifier* is not a valid path or model instance.

    Examples
    --------
    >>> resolve_hallucination("none", None) is None
    True
    >>> checker = resolve_hallucination("heuristic", None)  # doctest: +SKIP
    >>> callable(checker)                                    # doctest: +SKIP
    True
    >>> checker = resolve_hallucination("ml", "model_out/") # doctest: +SKIP
    """
    if mode == "none":
        return None
    if mode == "heuristic":
        from deepretro.algorithms.pipeline_checks import hallucination_checker
        return hallucination_checker
    if mode == "ml":
        from deepretro.models.hallucination_classifier import HallucinationClassifier

        if isinstance(classifier, (str, Path)):
            clf = HallucinationClassifier()
            clf.load(str(classifier))
        elif hasattr(classifier, "predict_single"):
            clf = classifier
        else:
            raise ValueError(
                f"hallucination_mode='ml' requires a HallucinationClassifier "
                f"or path to saved model — got {type(classifier)}"
            )
        return build_ml_checker(clf)
    raise ValueError(
        f"hallucination_mode must be 'heuristic', 'ml', or 'none' — got {mode!r}"
    )
