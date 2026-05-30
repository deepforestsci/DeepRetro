"""Helpers for optional hallucination filtering in retrosynthesis pipelines."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.typing import HallucinationChecker


def _as_pathway(candidate: Pathway | str) -> Pathway:
    """Normalize one candidate into a pathway list."""
    return [candidate] if isinstance(candidate, str) else list(candidate)


def build_ml_checker(classifier: Any) -> HallucinationChecker:
    """Wrap a classifier in the checker interface used by ``AutoSolver``.

    Parameters
    ----------
    classifier : Any
        Object exposing ``predict_single(product_smiles, reactants_smiles)``.

    Returns
    -------
    HallucinationChecker
        Callable returning ``(200, retained_pathways)``.
    """

    def _checker(product: str, pathways: list[Pathway]) -> tuple[int, list[Pathway]]:
        from deepretro.utils.utils_molecule import is_valid_smiles

        retained: list[Pathway] = []
        for candidate in pathways:
            pathway = _as_pathway(candidate)
            reactants_smiles = ".".join(pathway)
            if not is_valid_smiles(reactants_smiles):
                continue

            prediction = classifier.predict_single(product, reactants_smiles)
            if not bool(prediction.get("is_hallucination", True)):
                retained.append(pathway)

        return 200, retained

    return _checker


def filter_with_checker(
    product: str,
    pathways: list[Pathway],
    explanations: Iterable[str],
    confidence: Iterable[float],
    checker: HallucinationChecker | None,
) -> tuple[list[Pathway], list[str], list[float]]:
    """Filter pathways and keep explanations/confidence aligned.

    Parameters
    ----------
    product : str
        Product molecule SMILES.
    pathways : list[list[str]]
        Candidate precursor pathways.
    explanations : Iterable[str]
        Explanations aligned with ``pathways``.
    confidence : Iterable[float]
        Confidence scores aligned with ``pathways``.
    checker : HallucinationChecker or None
        Optional pathway filter.

    Returns
    -------
    tuple[list[list[str]], list[str], list[float]]
        Retained pathways with aligned metadata.
    """
    if checker is None or not pathways:
        return pathways, list(explanations), list(confidence)

    status, retained_pathways = checker(product, pathways)
    if status != 200:
        return [], [], []

    retained_explanations: list[str] = []
    retained_confidence: list[float] = []
    available = list(zip(pathways, explanations, confidence))
    for retained in retained_pathways:
        for index, (pathway, explanation, score) in enumerate(available):
            if pathway == retained:
                retained_explanations.append(explanation)
                retained_confidence.append(float(score))
                available.pop(index)
                break

    return retained_pathways, retained_explanations, retained_confidence


def resolve_hallucination(
    mode: str,
    classifier: str | Path | Any | None,
) -> HallucinationChecker | None:
    """Resolve an AutoSolver hallucination mode into a checker callable.

    Parameters
    ----------
    mode : str
        One of ``"heuristic"``, ``"ml"``, or ``"none"``.
    classifier : str, Path, object, or None
        Saved classifier directory or classifier object for ``"ml"`` mode.

    Returns
    -------
    HallucinationChecker or None
        ``None`` for ``"none"`` and ``"heuristic"`` because the built-in
        heuristic path is handled by ``deepretro.utils.llm.llm_pipeline``.
    """
    normalized_mode = mode.lower()
    if normalized_mode in {"none", "heuristic"}:
        return None
    if normalized_mode != "ml":
        raise ValueError(
            f"hallucination_mode must be 'heuristic', 'ml', or 'none' (got {mode!r})"
        )

    if isinstance(classifier, (str, Path)):
        from deepretro.models.hallucination_classifier import HallucinationClassifier

        loaded_classifier = HallucinationClassifier()
        loaded_classifier.load(str(classifier))
        return build_ml_checker(loaded_classifier)

    if hasattr(classifier, "predict_single"):
        return build_ml_checker(classifier)

    raise ValueError(
        "hallucination_mode='ml' requires a HallucinationClassifier instance "
        "or path to a saved classifier"
    )
