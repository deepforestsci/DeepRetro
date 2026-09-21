"""Optional adapter for official URSA route annotations.

URSA evaluates complete retrosynthetic routes. This module converts
DeepRetro's completed viewer-style route once, delegates every chemistry
decision to the official URSA/ChemCensor workflow, and returns the released
``PathResult.to_dict()`` payload unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem
from ursa import RetrosyntheticPath, Ursa

from deepretro._pathway import (
    _extract_step_smiles,
    _get_pathway_steps,
)


class UrsaError(RuntimeError):
    """Base error for optional URSA route annotation."""


class UrsaConfigurationError(UrsaError):
    """URSA asset configuration or backend initialization is invalid."""


class UrsaRouteError(UrsaError):
    """A DeepRetro pathway cannot be represented as one URSA route tree."""


class UrsaEvaluationError(UrsaError):
    """The official URSA backend failed while evaluating a converted route."""


@dataclass(frozen=True)
class UrsaConfig:
    """Configuration for optional URSA route annotation.

    Asset paths are required only when ``enabled`` is true. Requiring explicit
    paths prevents URSA from implicitly downloading its approximately 500 MB
    ChemCensor archive (about 5.7 GB extracted) during an ordinary DeepRetro
    run.
    """

    enabled: bool = False
    bb_catalog_path: Path | str | None = None
    chemcensor_db_path: Path | str | None = None

    def __post_init__(self) -> None:
        """Normalize and, when enabled, validate both official data assets."""
        catalog = _normalize_path(self.bb_catalog_path)
        database = _normalize_path(self.chemcensor_db_path)
        object.__setattr__(self, "bb_catalog_path", catalog)
        object.__setattr__(self, "chemcensor_db_path", database)

        if not self.enabled:
            return
        if catalog is None or not catalog.is_file():
            raise UrsaConfigurationError(
                "URSA building-block catalog was not found; set "
                "bb_catalog_path to the official URSA catalog file."
            )
        if database is None or not database.is_file():
            raise UrsaConfigurationError(
                "URSA ChemCensor database was not found; set "
                "chemcensor_db_path to the extracted official SQLite file."
            )


Pathway = dict[str, Any] | list[dict[str, Any]]


class UrsaAdapter:
    """Evaluate completed DeepRetro routes with official URSA."""

    def __init__(self, config: UrsaConfig) -> None:
        self.config = config
        self._backend: Ursa | None = None

    def annotate(
        self,
        pathway: Pathway,
        *,
        path_id: str = "deepretro-route",
    ) -> dict[str, Any]:
        """Return the official serialized URSA result for one completed route."""
        if not self.config.enabled:
            raise UrsaConfigurationError("URSA annotation is disabled.")

        route_data = _to_retrocast_route(pathway)
        backend = self._load_backend()

        try:
            official_path = RetrosyntheticPath.from_dict(route_data, path_id)
        except Exception as exc:
            raise UrsaRouteError(
                f"Official URSA route construction failed for {path_id!r}: {exc}"
            ) from exc

        try:
            payload = backend.score(official_path).to_dict()
        except Exception as exc:
            raise UrsaEvaluationError(
                f"URSA evaluation failed for {path_id!r}: {exc}"
            ) from exc
        if not isinstance(payload, dict):
            raise UrsaEvaluationError(
                f"URSA evaluation failed for {path_id!r}: "
                "PathResult.to_dict() did not return a dictionary."
            )
        return payload

    def _load_backend(self) -> Ursa:
        """Construct the official backend once."""
        if self._backend is None:
            try:
                self._backend = Ursa(
                    bb_catalog_path=self.config.bb_catalog_path,
                    chemcensor_db_path=self.config.chemcensor_db_path,
                    parallel=False,
                )
            except Exception as exc:
                raise UrsaConfigurationError(
                    f"Failed to initialize URSA with the configured assets: {exc}"
                ) from exc
        return self._backend


def _normalize_path(value: Path | str | None) -> Path | None:
    if value is None:
        return None
    return Path(value).expanduser().resolve()


def _to_retrocast_route(pathway: Pathway) -> dict[str, Any]:
    """Convert one flat DeepRetro route to URSA's recursive input shape."""
    steps, pathway_error = _get_pathway_steps(pathway)
    if pathway_error is not None:
        raise UrsaRouteError(pathway_error)
    if not steps:
        raise UrsaRouteError("URSA annotation requires at least one step.")

    producers: dict[str, tuple[str, list[str]]] = {}
    for index, step in enumerate(steps):
        if isinstance(step, dict) and "products" in step:
            products = step["products"]
            if isinstance(products, list) and len(products) != 1:
                raise UrsaRouteError(
                    f"Malformed step {index}: URSA requires exactly one product."
                )
        elif isinstance(step, dict):
            if "product" in step and "product_smiles" in step:
                raise UrsaRouteError(
                    f"Malformed step {index}: ambiguous product fields."
                )
            if "reactants" in step and "reactants_smiles" in step:
                raise UrsaRouteError(
                    f"Malformed step {index}: ambiguous reactant fields."
                )
        try:
            product, reactants = _extract_step_smiles(step)
        except ValueError as exc:
            raise UrsaRouteError(f"Malformed step {index}: {exc}") from exc

        product_key = _canonicalize_smiles(product)
        if product_key is None:
            raise UrsaRouteError(f"Invalid product SMILES at step {index}: {product!r}")
        normalized_reactants: list[str] = []
        for reactant in reactants:
            if _canonicalize_smiles(reactant) is None:
                raise UrsaRouteError(
                    f"Invalid reactant SMILES at step {index}: {reactant!r}"
                )
            normalized_reactants.append(reactant)

        if product_key in producers:
            raise UrsaRouteError(
                "Duplicate producing steps for canonical product "
                f"{product_key!r}; URSA requires an unambiguous route tree."
            )
        producers[product_key] = (product, normalized_reactants)

    consumed_products: set[str] = set()
    for _, reactants in producers.values():
        for reactant in reactants:
            reactant_key = _canonicalize_smiles(reactant)
            if reactant_key not in producers:
                continue
            if reactant_key in consumed_products:
                raise UrsaRouteError(
                    f"Produced intermediate {reactant_key!r} is consumed more than "
                    "once; URSA requires a route tree."
                )
            consumed_products.add(reactant_key)
    roots = [
        product_key for product_key in producers if product_key not in consumed_products
    ]
    if not roots:
        raise UrsaRouteError("URSA route contains a cycle or has no root product.")
    if len(roots) != 1:
        raise UrsaRouteError(
            f"URSA annotation requires exactly one root product; found {len(roots)}."
        )

    active: set[str] = set()
    visited: set[str] = set()

    def build_node(product_key: str) -> dict[str, Any]:
        if product_key in active:
            raise UrsaRouteError(
                f"URSA route contains a cycle involving {product_key!r}."
            )
        active.add(product_key)
        visited.add(product_key)
        product, reactants = producers[product_key]
        child_nodes: list[dict[str, Any]] = []
        for reactant in reactants:
            reactant_key = _canonicalize_smiles(reactant)
            if reactant_key in producers:
                child_nodes.append(build_node(reactant_key))
            else:
                child_nodes.append({"smiles": reactant, "product_of": None})
        active.remove(product_key)
        return {
            "smiles": product,
            "product_of": {"reactants": child_nodes},
        }

    target = build_node(roots[0])
    if visited != producers.keys():
        raise UrsaRouteError(
            "URSA route contains disconnected or cyclic producing steps."
        )
    return {"schema_version": "2", "target": target}


def _canonicalize_smiles(smiles: str) -> str | None:
    """Return canonical isomeric SMILES for URSA route graph linkage."""
    text = str(smiles).strip() if smiles is not None else ""
    if not text:
        return None
    try:
        molecule = Chem.MolFromSmiles(text)
    except Exception:
        return None
    if molecule is None:
        return None
    return str(Chem.MolToSmiles(molecule, isomericSmiles=True))


__all__ = [
    "UrsaAdapter",
    "UrsaConfig",
    "UrsaConfigurationError",
    "UrsaError",
    "UrsaEvaluationError",
    "UrsaRouteError",
]
