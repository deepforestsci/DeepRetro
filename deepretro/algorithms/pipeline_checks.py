"""Pipeline-compatible wrappers for package-native heuristics."""

from __future__ import annotations

import structlog

from deepretro.algorithms.hallucination_checker import calculate_hallucination_score
from deepretro.algorithms.hallucination_weights import (
    DEFAULT_WEIGHTS,
    HallucinationWeights,
)
from deepretro.utils.utils_molecule import is_valid_smiles

logger = structlog.get_logger(__name__)

ALLOWED_HALLUCINATION_SEVERITIES = {"low", "medium"}


def hallucination_checker(
    product: str,
    res_smiles: list[str | list[str]] | str,
    weights: HallucinationWeights | None = None,
    rank_only: bool = True,
) -> tuple[int, list[list[str]]]:
    """Rank candidate pathways by hallucination score, least suspicious first.

    Assessable candidates are ordered by descending score, then ascending raw
    penalty total. Flagged candidates remain fallback options by default.
    Retention is not an individual chemical-validity verdict. Invalid, empty,
    no-op and carbon-radical inputs are unsupported and excluded in both modes.

    Parameters
    ----------
    product : str
        SMILES string of the target product molecule.
    res_smiles : list
        Candidate pathways represented as either lists of reactant SMILES
        strings or a single reactant SMILES string. A bare string is treated as
        one pathway, not iterated character by character.
    weights : HallucinationWeights, optional
        Penalty weights and cutoffs. ``None`` (the default) uses the default
        configuration. The reject_below threshold is independent of severity.
    rank_only : bool, optional
        When ``True`` (the default), return every assessable pathway,
        accepted ones first, each tier ordered by descending score. When
        ``False``, use the historical hard-filter behaviour and return
        ``(400, [])`` if nothing is accepted.

    Returns
    -------
    tuple[int, list]
        ``(200, pathways)`` with accepted pathways first. Under ``rank_only``
        the status is 400 when no assessable pathway remains.

    Examples
    --------
    >>> status, ranked = hallucination_checker("CC=O", [["CCO"]])
    >>> status
    200
    >>> ranked
    [['CCO']]
    """
    # A bare SMILES string is a single pathway. Without this it is iterated
    # character by character and each character becomes its own "pathway".
    if isinstance(res_smiles, str):
        res_smiles = [res_smiles]

    accepted: list[tuple[int, int, list[str]]] = []
    rejected: list[tuple[int, int, list[str], str]] = []
    unassessable: list[tuple[list[str], str]] = []

    for pathway in res_smiles:
        if isinstance(pathway, list):
            reactants_smi = ".".join(pathway)
            normalized_pathway = pathway
        else:
            reactants_smi = pathway
            normalized_pathway = [pathway]

        if not is_valid_smiles(reactants_smi):
            continue

        report = calculate_hallucination_score(reactants_smi, product, weights)
        # Gate on the score, not the severity label. The label collapses four
        # levels into one boundary anyway, and going through it tied the gate to
        # the cut_high < cut_medium invariant.
        threshold = (weights or DEFAULT_WEIGHTS).reject_below

        if report.get("unassessable"):
            unassessable.append((normalized_pathway, report["message"]))
            continue

        # penalty_total orders the demoted tier below zero, where `score` floors
        # and 51% of demoted candidates tie.
        penalty = report.get("penalty_total", 0)
        if report["score"] >= threshold:
            accepted.append((report["score"], penalty, normalized_pathway))
        else:
            rejected.append(
                (report["score"], penalty, normalized_pathway, report["severity"])
            )

    if not rank_only:
        if accepted:
            return 200, [pathway for _, _, pathway in accepted]
        return 400, []

    # Least suspicious first: descending score, then ascending penalty total to
    # order the candidates whose score has floored at 0.
    accepted.sort(key=lambda item: (-item[0], item[1]))
    rejected.sort(key=lambda item: (-item[0], item[1]))

    if unassessable:
        logger.info(
            "Hallucination check dropped unassessable candidates",
            product=product,
            n_dropped=len(unassessable),
            dropped=[
                {"reactants": ".".join(pathway), "reason": reason}
                for pathway, reason in unassessable
            ],
        )

    if rejected:
        # The verdict used to be computed and thrown away, so a chemist could
        # never tell that a route leaned on a flagged step.
        logger.info(
            "Hallucination check demoted candidates",
            product=product,
            n_accepted=len(accepted),
            n_demoted=len(rejected),
            demoted=[
                {
                    "reactants": ".".join(pathway),
                    "score": score,
                    "severity": severity,
                }
                for score, _penalty, pathway, severity in rejected
            ],
        )

    ranked = [pathway for _, _, pathway in accepted]
    ranked += [pathway for _, _, pathway, _ in rejected]
    if ranked:
        return 200, ranked
    return 400, []
