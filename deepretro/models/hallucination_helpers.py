"""Helpers for wiring hallucination checkers into the retrosynthesis pipeline.

Provides:

* :class:`MLChecker` — wrap a classifier into the callable
  signature the pipeline expects.
* :func:`resolve_hallucination` — turn a string to determine checker type
  into a checker callable (or ``None``) consumed by ``llm_pipeline``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

from deepretro.models.hallucination_classifier import predict_single_reaction
from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.typing import HallucinationChecker


class MLChecker:
    """Wrap a classifier into a pipeline-compatible hallucination checker.

    The checker has the same signature as the built-in heuristic checker
    (``(product, pathways) -> (int, list)``). Pathways flagged as
    hallucinated are dropped; if all are rejected the pipeline retries the
    LLM call.

    Parameters
    ----------
    clf : HallucinationClassifier
        A fitted classifier instance or classifier-like object with
        ``predict_probability()``, ``threshold``, and ``featurizer``
        attributes compatible with :func:`predict_single_reaction`.

    Examples
    --------
    >>> from deepretro.models import HallucinationClassifier  # doctest: +SKIP
    >>> clf = HallucinationClassifier()                       # doctest: +SKIP
    >>> clf.load("model_out/")                               # doctest: +SKIP
    >>> checker = MLChecker(clf)                             # doctest: +SKIP
    >>> status, kept = checker("CCO", [["CC", "O"]])         # doctest: +SKIP
    """

    def __init__(self, clf: Any) -> None:
        """Store the classifier and SMILES validator used by the checker.

        Parameters
        ----------
        clf : Any
            Classifier-like object implementing
            ``predict_probability()``, ``threshold``, and ``featurizer``
            attributes compatible with :func:`predict_single_reaction`.
        """
        from deepretro.utils.utils_molecule import is_valid_smiles

        self.clf = clf
        self.is_valid_smiles = is_valid_smiles

    def __call__(self, product: str, pathways: list) -> tuple[int, list]:
        """Filter pathway candidates using the configured ML classifier.

        Parameters
        ----------
        product : str
            SMILES string of the target product.
        pathways : list
            Candidate pathways represented as either lists of reactant
            SMILES strings or single dot-separated reactant strings.

        Returns
        -------
        tuple[int, list]
            ``(200, valid_pathways)`` when at least one pathway is kept,
            otherwise ``(400, [])``.
        """
        valid_pathways = []
        for pathway in pathways:
            if isinstance(pathway, list):
                reactants_smi = ".".join(pathway)
                normalized_pathway = pathway
            else:
                reactants_smi = pathway
                normalized_pathway = [pathway]

            if not self.is_valid_smiles(reactants_smi):
                continue

            pred = predict_single_reaction(self.clf, product, reactants_smi)
            if pred.get("is_hallucination", True):
                continue

            valid_pathways.append(normalized_pathway)

        if valid_pathways:
            return 200, valid_pathways
        return 400, []


def resolve_hallucination(
    mode: str,
    classifier: str | Path | Any | None,
) -> Callable | None:
    """Convert a user-facing hallucination mode into a checker callable.

    Parameters
    ----------
    mode : str
        One of ``"heuristic"``, ``"ml"``, ``"none"``, or ``"default"``.
        ``"default"`` returns the placeholder checker from
        :func:`deepretro.models.hallucination_checker.make_hallucination_checker`.
    classifier : str or Path or Any or None
        Path to a saved ML model directory, or a classifier-like object
        exposing ``predict_probability()``, ``threshold``, and
        ``featurizer`` attributes compatible with
        :func:`predict_single_reaction`. Required
        when *mode* is ``"ml"``, ignored otherwise.

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
    >>> from pathlib import Path
    >>> resolve_hallucination("none", None) is None
    True
    >>> checker = resolve_hallucination("heuristic", None)  # doctest: +SKIP
    >>> callable(checker)                                    # doctest: +SKIP
    True
    >>> from deepretro.models import HallucinationClassifier  # doctest: +SKIP
    >>> clf = HallucinationClassifier()                       # doctest: +SKIP
    >>> checker = resolve_hallucination("ml", clf)           # doctest: +SKIP
    >>> checker = resolve_hallucination("ml", Path("model_out/"))  # doctest: +SKIP
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
        elif hasattr(classifier, "predict_probability"):
            clf = classifier
        else:
            raise ValueError(
                f"hallucination_mode='ml' requires a HallucinationClassifier "
                f"or path to saved model — got {type(classifier)}"
            )
        return MLChecker(clf)
    if mode == "default":
        from deepretro.models.hallucination_checker import make_hallucination_checker

        return make_hallucination_checker()
    raise ValueError(
        f"hallucination_mode must be 'heuristic', 'ml', 'none', or 'default' "
        f"— got {mode!r}"
    )


def _as_pathway(candidate: Pathway | str) -> Pathway:
    """Normalize one candidate into a list-of-SMILES pathway.

    Parameters
    ----------
    candidate : list[str] or str
        A pathway as a list of reactant SMILES, or a single reactant SMILES.

    Returns
    -------
    list[str]
        The candidate as a list, so string and list pathways compare equal.
    """
    return [candidate] if isinstance(candidate, str) else list(candidate)


def filter_with_checker(
    product: str,
    pathways: list[Pathway],
    explanations: Iterable[str],
    confidence: Iterable[float],
    checker: HallucinationChecker | None,
) -> tuple[list[Pathway], list[str], list[float]]:
    """Filter pathways with a checker while keeping metadata aligned.

    A checker returns only the kept pathways, not their indices, so the
    explanations and confidence scores must be realigned to the surviving
    pathways. Alignment is by value (order-preserving), which is why a
    checker must only ever return pathways from its input.

    Parameters
    ----------
    product : str
        Product molecule SMILES.
    pathways : list[list[str]]
        Candidate precursor pathways.
    explanations : Iterable[str]
        Explanations aligned position-for-position with ``pathways``.
    confidence : Iterable[float]
        Confidence scores aligned position-for-position with ``pathways``.
    checker : HallucinationChecker or None
        Optional pathway filter ``(product, pathways) -> (status, kept)``.
        ``None`` disables filtering and returns the inputs unchanged.

    Returns
    -------
    tuple[list[list[str]], list[str], list[float]]
        Retained pathways with explanations and confidence realigned.
        A checker status other than 200 yields three empty lists.

    Raises
    ------
    ValueError
        If ``pathways``, ``explanations``, and ``confidence`` differ in
        length, or if the checker returns a pathway absent from the input.

    Examples
    --------
    >>> filter_with_checker(
    ...     "CC(=O)Oc1ccccc1C(=O)O",
    ...     [["OC(=O)c1ccccc1O"], ["CC(=O)O"]],
    ...     ["hydrolysis", "deacetylation"],
    ...     [0.3, 0.9],
    ...     None,
    ... )
    ([['OC(=O)c1ccccc1O'], ['CC(=O)O']], ['hydrolysis', 'deacetylation'], [0.3, 0.9])
    """
    explanation_list = list(explanations)
    confidence_list = list(confidence)
    if len(pathways) != len(explanation_list) or len(pathways) != len(confidence_list):
        raise ValueError(
            f"pathways ({len(pathways)}), explanations ({len(explanation_list)}), "
            f"and confidence ({len(confidence_list)}) must have the same length"
        )

    if checker is None or not pathways:
        return list(pathways), explanation_list, confidence_list

    status, retained_pathways = checker(product, pathways)
    if status != 200:
        return [], [], []

    available = [
        (_as_pathway(pathway), explanation, score)
        for pathway, explanation, score in zip(
            pathways, explanation_list, confidence_list
        )
    ]
    matched_pathways: list[Pathway] = []
    matched_explanations: list[str] = []
    matched_confidence: list[float] = []

    for retained in retained_pathways:
        target = _as_pathway(retained)
        for index, (normalized, explanation, score) in enumerate(available):
            if normalized == target:
                matched_pathways.append(retained)
                matched_explanations.append(explanation)
                matched_confidence.append(float(score))
                available.pop(index)
                break
        else:
            raise ValueError(
                f"Checker returned pathway {retained!r} not present in original input"
            )

    return matched_pathways, matched_explanations, matched_confidence
