"""Placeholder hallucination checker that pins the runtime callable contract.

This module reserves the model-backed slot in the retrosynthesis pipeline. It
is deliberately distinct from ``deepretro.algorithms.hallucination_checker``,
which is the heuristic (atom/ring/aromaticity) checker; per the repo's
``algorithms`` (heuristic) vs ``models`` (learned) split, a trained classifier
belongs here.

Nothing is wired into a live filtering path yet. The point of this file is to
lock in the callable seam — ``deepretro.utils.typing.HallucinationChecker``,
i.e. ``(product, pathways) -> (status, kept)`` — so a real model can replace
the body later without any caller renaming anything.

The verdict below carries no chemical meaning. It is a seeded hash bit, chosen
only so the callable is deterministic and both return shapes are exercised in
tests. Do not read anything into which molecules keep or flag.
"""

from __future__ import annotations

import hashlib

import structlog

from deepretro.utils.typing import HallucinationChecker

logger = structlog.get_logger(__name__)


def make_hallucination_checker(seed: int = 0) -> HallucinationChecker:
    """Build a deterministic placeholder checker.

    The returned callable implements the production contract
    ``(product, pathways) -> (status, kept)``:

    * keep -> ``(200, pathways)`` (pathways returned unchanged)
    * flag -> ``(400, [])`` (pathways rejected)

    A trained model can replace the inner body behind this same signature and
    return type; callers and the agent-loop tool wiring stay unchanged.
    """

    def check_hallucination(product: str, pathways: list) -> tuple[int, list]:
        # Placeholder verdict: a seeded hash bit, NOT a chemical judgement.
        digest = hashlib.sha256(f"{seed}|{product}|{pathways}".encode()).hexdigest()
        keep = int(digest, 16) % 2 == 1

        logger.info("hallucination_checker.verdict", product=product, keep=keep)

        if keep:
            return 200, pathways
        return 400, []

    return check_hallucination
