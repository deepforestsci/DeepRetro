"""Shared type aliases for the deepretro.utils package."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

RouteNode = Mapping[str, Any]
Step = dict[str, Any]
DependencyMap = dict[str, list[str]]
ParseOutput = dict[str, Any]

# A hallucination checker filters candidate pathways for one product.
# ``(product_smiles, pathways) -> (status_code, kept_pathways)`` where a
# status of 200 means at least one pathway was kept.
HallucinationChecker = Callable[[str, list], tuple[int, list]]
