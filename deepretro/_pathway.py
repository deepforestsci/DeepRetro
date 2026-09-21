"""Shared parsing helpers for DeepRetro pathway step shapes."""

from __future__ import annotations

from typing import Any


def _get_pathway_steps(pathway: Any) -> tuple[list[Any], str | None]:
    if isinstance(pathway, dict):
        if "steps" not in pathway:
            return [], "Malformed pathway: missing steps"
        steps = pathway["steps"]
        if not isinstance(steps, list):
            return [], "Malformed pathway: steps must be a list"
        return steps, None
    if isinstance(pathway, list):
        return pathway, None
    return [], "Malformed pathway: expected dict or list"


def _extract_step_smiles(step: Any) -> tuple[str, list[str]]:
    if not isinstance(step, dict):
        raise ValueError("step must be a dictionary")

    if "products" in step:
        if "reactants" not in step:
            raise ValueError("mixed or incomplete viewer-style step")
        if "product" in step or "product_smiles" in step or "reactants_smiles" in step:
            raise ValueError("mixed viewer-style and simple-style step")
        return _extract_viewer_step_smiles(step)

    has_simple_product = "product_smiles" in step or "product" in step
    has_simple_reactants = "reactants_smiles" in step or "reactants" in step
    if has_simple_product and has_simple_reactants:
        product = step.get("product_smiles", step.get("product"))
        reactants = step.get("reactants_smiles", step.get("reactants"))
        if product is None:
            raise ValueError("missing product SMILES")
        reactant_list = _split_reactants(reactants)
        if not reactant_list:
            raise ValueError("missing reactant SMILES")
        return str(product), reactant_list

    raise ValueError("missing product or reactant SMILES")


def _extract_viewer_step_smiles(step: dict[str, Any]) -> tuple[str, list[str]]:
    products = step["products"]
    reactants = step["reactants"]
    if not isinstance(products, list) or not products:
        raise ValueError("missing product SMILES")
    if not isinstance(reactants, list):
        raise ValueError("missing reactant SMILES")

    product = _smiles_from_item(products[0])
    reactant_list = [_smiles_from_item(item).strip() for item in reactants]
    reactant_list = [smiles for smiles in reactant_list if smiles]
    if not reactant_list:
        raise ValueError("missing reactant SMILES")
    return product, reactant_list


def _split_reactants(reactants_smiles: Any) -> list[str]:
    if isinstance(reactants_smiles, str):
        parts = reactants_smiles.split(".")
    elif isinstance(reactants_smiles, list):
        parts = [str(item) for item in reactants_smiles]
    else:
        return []
    return [part.strip() for part in parts if part.strip()]


def _smiles_from_item(item: Any) -> str:
    value = item.get("smiles") if isinstance(item, dict) else item
    if value is None or not str(value).strip():
        raise ValueError("missing SMILES")
    return str(value)
