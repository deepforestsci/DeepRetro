"""Helpers for wiring hallucination checkers into the retrosynthesis pipeline.

Provides:

* :func:`build_ml_checker` — wrap a classifier into the callable
  signature the pipeline expects.
* :func:`resolve_hallucination_args` — turn a user-friendly mode string
  into the ``(hallucination_check, hallucination_checker_fn)`` pair
  consumed by ``llm_pipeline``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from deepretro.utils.utils_molecule import is_valid_smiles

VALID_MODES = ("heuristic", "ml", "none")


def build_ml_checker(clf: Any):
    """Wrap a ``HallucinationClassifier`` into a callable with the same
    signature as ``src.utils.hallucination_checks.hallucination_checker``:

        (product: str, pathways: list) -> (int, list)

    Pathways flagged as hallucinated are dropped, exactly like the
    heuristic checker.  This plugs into ``llm_pipeline``'s retry loop
    so rejected results trigger a new LLM call.
    """
    def _checker(product: str, pathways: list) -> tuple[int, list]:
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
) -> tuple[str, Any]:
    """Return ``(hallucination_check, hallucination_checker_fn)`` for the
    pipeline based on *hallucination_mode*.

    Parameters
    ----------
    hallucination_mode : str
        One of ``"heuristic"``, ``"ml"``, or ``"none"``.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when *hallucination_mode* is ``"ml"``.  Pass a fitted
        ``HallucinationClassifier`` instance or a ``str`` / ``Path``
        pointing to a saved model directory.
    """
    if hallucination_mode not in VALID_MODES:
        raise ValueError(
            f"hallucination_mode must be one of {VALID_MODES}, "
            f"got {hallucination_mode!r}"
        )

    if hallucination_mode == "none":
        return "False", None

    if hallucination_mode == "heuristic":
        return "True", None

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

    return "True", build_ml_checker(clf)
