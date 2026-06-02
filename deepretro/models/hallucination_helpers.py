"""Helpers for wiring hallucination checkers into the retrosynthesis pipeline.

Provides:

* :class:`MLChecker` — wrap a classifier into the callable
  signature the pipeline expects.
* :func:`resolve_hallucination` — turn a string to determine checker type
  into a checker callable (or ``None``) consumed by ``llm_pipeline``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from deepretro.models.hallucination_classifier import predict_single_reaction


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
        One of ``"heuristic"``, ``"ml"``, or ``"none"``.
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
    raise ValueError(
        f"hallucination_mode must be 'heuristic', 'ml', or 'none' — got {mode!r}"
    )
