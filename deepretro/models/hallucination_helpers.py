"""Helpers for wiring hallucination checkers into the retrosynthesis pipeline.

Provides:

* :func:`build_ml_checker` — wrap a classifier into the callable
  signature the pipeline expects.
* :func:`resolve_hallucination_args` — turn a user-friendly mode string
  into a checker callable (or ``None``) consumed by ``llm_pipeline``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

VALID_MODES = ("heuristic", "ml", "none")


def build_ml_checker(clf: Any) -> Any:
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
    checker : callable
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


def resolve_hallucination_args(
    hallucination_mode: str,
    hallucination_classifier: Any,
):
    """Translate a user-facing mode string into a checker callable or None.

    Parameters
    ----------
    hallucination_mode : str
        One of ``"heuristic"``, ``"ml"``, or ``"none"``.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when *hallucination_mode* is ``"ml"``.  Pass a fitted
        ``HallucinationClassifier`` instance or a ``str`` / ``Path``
        pointing to a saved model directory.

    Returns
    -------
    callable or None
        A checker with signature ``(product, pathways) -> (int, list)``,
        or ``None`` to skip hallucination checking.

    Raises
    ------
    ValueError
        If *hallucination_mode* is not one of the valid modes, or if
        ``"ml"`` mode is requested without a valid classifier.

    Examples
    --------
    >>> resolve_hallucination_args("none", None) is None
    True
    >>> callable(resolve_hallucination_args("heuristic", None))
    True
    """
    if hallucination_mode not in VALID_MODES:
        raise ValueError(
            f"hallucination_mode must be one of {VALID_MODES}, "
            f"got {hallucination_mode!r}"
        )

    if hallucination_mode == "none":
        return None

    if hallucination_mode == "heuristic":
        from deepretro.algorithms.pipeline_checks import hallucination_checker
        return hallucination_checker

    # mode == "ml" -- resolve the classifier
    from deepretro.models.hallucination_classifier import HallucinationClassifier

    if isinstance(hallucination_classifier, (str, Path)):
        clf = HallucinationClassifier()
        clf.load(str(hallucination_classifier))
    elif isinstance(hallucination_classifier, HallucinationClassifier):
        clf = hallucination_classifier
    else:
        raise ValueError(
            "hallucination_mode='ml' requires hallucination_classifier "
            "to be a HallucinationClassifier instance or a path to a "
            f"saved model directory — got {type(hallucination_classifier)}"
        )

    return build_ml_checker(clf)
