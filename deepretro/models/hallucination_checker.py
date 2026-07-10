"""Temporary deterministic hallucination-checker skeleton."""

from __future__ import annotations

import hashlib

import structlog

from deepretro.utils.typing import HallucinationChecker

logger = structlog.get_logger(__name__)


def make_hallucination_checker(seed: int = 0) -> HallucinationChecker:
    """Build a deterministic checker callable."""

    def check_hallucination(product: str, pathways: list) -> tuple[int, list]:
        """Return a deterministic keep or flag verdict."""
        digest = hashlib.sha256(f"{seed}|{product}|{pathways}".encode()).hexdigest()
        verdict = int(digest, 16) % 2
        keep = verdict == 1

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
