"""Helpers for optional hallucination filtering in retrosynthesis pipelines."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Protocol, runtime_checkable

from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.typing import HallucinationChecker
from deepretro.utils.utils_molecule import is_valid_smiles


@runtime_checkable
class HallucinationPredictor(Protocol):
    """Protocol for hallucination classifiers used by :func:`build_ml_checker`."""

    def predict_single(
        self, product_smiles: str, reactants_smiles: str
    ) -> dict[str, Any]: ...


def _as_pathway(candidate: Pathway | str) -> Pathway:
    """Normalize one candidate into a pathway list."""
    return [candidate] if isinstance(candidate, str) else list(candidate)


def build_ml_checker(classifier: HallucinationPredictor) -> HallucinationChecker:
    """Wrap a classifier in the checker interface used by ``AutoSolver``.

    Parameters
    ----------
    classifier : HallucinationPredictor
        Object exposing ``predict_single(product_smiles, reactants_smiles)``.

    Returns
    -------
    HallucinationChecker
        Callable returning ``(200, retained_pathways)``.

    Raises
    ------
    TypeError
        If *classifier* does not expose ``predict_single``.

    Examples
    --------
    >>> class MyClassifier:
    ...     def predict_single(self, product, reactants):
    ...         return {"is_hallucination": False}
    >>> checker = build_ml_checker(MyClassifier())
    >>> status, retained = checker("CC(=O)Oc1ccccc1C(=O)O", [["OC(=O)c1ccccc1O"]])
    >>> status
    200
    >>> retained
    [['OC(=O)c1ccccc1O']]
    """
    if not callable(getattr(classifier, "predict_single", None)):
        raise TypeError(
            "classifier must expose predict_single(product_smiles, reactants_smiles)"
        )

    def _checker(product: str, pathways: list[Pathway]) -> tuple[int, list[Pathway]]:
        retained: list[Pathway] = []
        for candidate in pathways:
            pathway = _as_pathway(candidate)
            reactants_smiles = ".".join(pathway)
            if not is_valid_smiles(reactants_smiles):
                continue

            prediction = classifier.predict_single(product, reactants_smiles)
            if prediction.get("is_hallucination") is False:
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
        Retained pathways with aligned metadata. Only pathways present in
        the original input are included; unknown pathways returned by the
        checker raise ``ValueError``.

    Raises
    ------
    ValueError
        If ``pathways``, ``explanations``, and ``confidence`` have
        different lengths, or if the checker returns a pathway not present
        in the original input.

    Examples
    --------
    >>> pathways, expl, conf = filter_with_checker(
    ...     "CC(=O)Oc1ccccc1C(=O)O",
    ...     [["OC(=O)c1ccccc1O"], ["CC(=O)O"]],
    ...     ["hydrolysis", "deacetylation"],
    ...     [0.3, 0.9],
    ...     None,
    ... )
    >>> pathways
    [['OC(=O)c1ccccc1O'], ['CC(=O)O']]
    >>> conf
    [0.3, 0.9]
    """
    explanation_list = list(explanations)
    confidence_list = list(confidence)
    if len(pathways) != len(explanation_list) or len(pathways) != len(confidence_list):
        raise ValueError(
            f"pathways ({len(pathways)}), explanations ({len(explanation_list)}), "
            f"and confidence ({len(confidence_list)}) must have the same length"
        )

    if checker is None or not pathways:
        return pathways, explanation_list, confidence_list

    status, retained_pathways = checker(product, pathways)
    if status != 200:
        return [], [], []

    available = list(zip(pathways, explanation_list, confidence_list))
    matched_pathways: list[Pathway] = []
    matched_explanations: list[str] = []
    matched_confidence: list[float] = []

    for retained in retained_pathways:
        found = False
        for index, (pathway, explanation, score) in enumerate(available):
            if pathway == retained:
                matched_pathways.append(retained)
                matched_explanations.append(explanation)
                matched_confidence.append(float(score))
                available.pop(index)
                found = True
                break
        if not found:
            raise ValueError(
                f"Checker returned pathway {retained!r} not present in original input"
            )

    return matched_pathways, matched_explanations, matched_confidence


def resolve_hallucination(
    mode: str,
    classifier: HallucinationPredictor | None,
) -> HallucinationChecker | None:
    """Resolve an AutoSolver hallucination mode into a checker callable.

    Parameters
    ----------
    mode : str
        One of ``"heuristic"``, ``"ml"``, or ``"none"``.
    classifier : HallucinationPredictor or None
        Classifier object exposing ``predict_single`` for ``"ml"`` mode.
        The caller is responsible for loading the model before passing it.

    Returns
    -------
    HallucinationChecker or None
        ``None`` for ``"none"`` and ``"heuristic"`` because the built-in
        heuristic path is handled by ``deepretro.utils.llm.llm_pipeline``.

    Raises
    ------
    ValueError
        If *mode* is not one of the three valid options, or if ``"ml"``
        mode is requested without a valid classifier.

    Examples
    --------
    >>> resolve_hallucination("none", None) is None
    True
    >>> resolve_hallucination("heuristic", None) is None
    True
    >>> class MyClassifier:
    ...     def predict_single(self, product, reactants):
    ...         return {"is_hallucination": False}
    >>> checker = resolve_hallucination("ml", MyClassifier())
    >>> callable(checker)
    True
    """
    normalized_mode = mode.lower()
    if normalized_mode in {"none", "heuristic"}:
        return None
    if normalized_mode != "ml":
        raise ValueError(
            f"hallucination_mode must be 'heuristic', 'ml', or 'none' (got {mode!r})"
        )

    if isinstance(classifier, HallucinationPredictor):
        return build_ml_checker(classifier)

    raise ValueError(
        "hallucination_mode='ml' requires a classifier object exposing "
        "predict_single(product_smiles, reactants_smiles)"
    )
