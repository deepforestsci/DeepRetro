"""Deterministic hallucination checker for the retrosynthesis agent.

:func:`make_hallucination_checker` builds a callable matching
:data:`~deepretro.utils.typing.HallucinationChecker` — ``(product, pathways) ->
(status, kept)`` — that the agent loop exposes as the ``check_hallucination``
tool. The verdict is a stable function of the inputs, so the same
``(product, pathways, seed)`` always produces the same decision.

The current body derives its keep/flag verdict from a seeded hash of the inputs
(a temporary deterministic placeholder to be replaced by a trained model).
Because every public name here describes the *role* — hallucination checking —
and never the mechanism, the internals can later be swapped for a real model
without any caller having to rename anything.
"""

from __future__ import annotations

import hashlib

import structlog

from deepretro.utils.tracing import observe
from deepretro.utils.typing import HallucinationChecker

logger = structlog.get_logger(__name__)


def make_hallucination_checker(seed: int = 0) -> HallucinationChecker:
    """Build a deterministic hallucination checker callable.

    Parameters
    ----------
    seed : int, optional
        Seed folded into the verdict hash so distinct seeds yield distinct but
        reproducible decision streams. Defaults to ``0``.

    Returns
    -------
    HallucinationChecker
        A callable ``(product: str, pathways: list) -> tuple[int, list]`` that
        returns ``(200, pathways)`` to keep the candidates or ``(400, [])`` to
        flag them as a hallucination. The verdict is deterministic in
        ``(product, pathways, seed)``.

    Notes
    -----
    The verdict is a temporary deterministic placeholder: it hashes the inputs
    to a value in ``{0, 1}`` (``1`` keeps, ``0`` flags). A trained classifier is
    expected to replace this body without changing the returned callable's
    contract.

    Examples
    --------
    >>> checker = make_hallucination_checker(seed=0)
    >>> keep = checker("CC=O", [["CCO"]])  # doctest: +ELLIPSIS
    [hallucination-checker] check_hallucination called: product='CC=O' reactants=[['CCO']] -> keep=True...
    >>> keep
    (200, [['CCO']])
    >>> flag = checker("CCO", [["CC", "O"]])  # doctest: +ELLIPSIS
    [hallucination-checker] check_hallucination called: product='CCO' reactants=[['CC', 'O']] -> keep=False...
    >>> flag
    (400, [])
    """

    @observe(name="check_hallucination")
    def check_hallucination(product: str, pathways: list) -> tuple[int, list]:
        """Return a deterministic keep/flag verdict for a proposed step.

        Parameters
        ----------
        product : str
            Product molecule SMILES.
        pathways : list
            Candidate precursor pathways for the product.

        Returns
        -------
        tuple[int, list]
            ``(200, pathways)`` to keep the candidates, ``(400, [])`` to flag
            them as a hallucination.
        """
        # Temporary deterministic placeholder for a trained model: fold the
        # inputs into a stable hash and read off a single bit as the verdict.
        digest = hashlib.sha256(f"{seed}|{product}|{pathways}".encode()).hexdigest()
        verdict = int(digest, 16) % 2
        keep = verdict == 1

        print(
            f"[hallucination-checker] check_hallucination called: "
            f"product={product!r} reactants={pathways} -> keep={keep}"
        )
        logger.info(
            "hallucination_checker.verdict",
            product=product,
            pathways=pathways,
            keep=keep,
        )

        if keep:
            return 200, pathways
        return 400, []

    return check_hallucination
