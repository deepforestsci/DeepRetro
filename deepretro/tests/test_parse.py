"""Unit tests for route tree parsing utilities."""

from __future__ import annotations

import inspect

from deepretro.utils.parse import RetrosynthesisRouteParser


def test_private_route_parser_helpers_are_documented_and_typed() -> None:
    """Private parser helpers should still have docstrings and type annotations."""
    private_helpers = [
        member
        for name, member in inspect.getmembers(
            RetrosynthesisRouteParser, inspect.isfunction
        )
        if name.startswith("_") and not name.startswith("__")
    ]

    assert private_helpers
    for helper in private_helpers:
        signature = inspect.signature(helper)
        assert inspect.getdoc(helper), helper.__name__
        assert all(
            parameter.annotation is not inspect.Signature.empty
            for parameter in signature.parameters.values()
            if parameter.name not in {"self", "cls"}
        ), helper.__name__
        assert signature.return_annotation is not inspect.Signature.empty, (
            helper.__name__
        )


def test_route_parser_formats_single_reaction_step_without_confidence() -> None:
    """Route parser should format steps without confidence estimates."""
    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=lambda smiles: f"formula:{smiles}",
        mass_calculator=lambda smiles: float(len(smiles)),
        scalability_calculator=lambda reactant, product: f"{reactant}->{product}",
    )
    route = {
        "smiles": "CCO",
        "children": [
            {
                "children": [{"smiles": "CC", "is_reaction": False}],
            }
        ],
    }

    result = parser.format_output(route)

    assert result["dependencies"] == {"1": []}
    assert result["steps"][0]["reactionmetrics"] == [
        {
            "scalabilityindex": "CC->CCO",
            "closestliterature": "",
        }
    ]


def test_route_parser_rebuilds_dependencies_from_step_products() -> None:
    """Route parser should derive step dependencies from product/reactant matches."""
    parser = RetrosynthesisRouteParser()
    step_list = [
        {"step": "1", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "B"}], "reactants": [{"smiles": "A"}]},
        {"step": "3", "products": [{"smiles": "C"}], "reactants": [{"smiles": "B"}]},
    ]

    assert parser.fix_dependencies(step_list) == {
        "1": [],
        "2": ["1"],
        "3": ["2"],
    }
