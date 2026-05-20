from collections.abc import Callable

import pytest

from deepretro.utils.parse import (
    RetrosynthesisRouteParser,
    fix_dependencies,
    format_output,
    parse_step,
)


def test_format_output_parses_single_step_with_metadata(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """Golden path: a one-step route (ethanol from ethane + water) should
    produce the correct step structure, dependency map, and real chemical
    metadata (formula C2H6O, mass ~46.04) computed by RDKit."""
    data = {
        "smiles": "CCO",
        "children": [{"children": [{"smiles": "CC"}, {"smiles": "O"}]}],
    }

    output = route_parser_factory().format_output(data)

    assert output["dependencies"] == {"1": []}
    assert len(output["steps"]) == 1

    step = output["steps"][0]
    assert step["step"] == "1"
    assert [p["smiles"] for p in step["products"]] == ["CCO"]
    assert [r["smiles"] for r in step["reactants"]] == ["CC", "O"]
    assert step["reactionmetrics"][0]["scalabilityindex"] == 0
    assert step["products"][0]["product_metadata"]["chemical_formula"] == "C2H6O"
    assert step["products"][0]["product_metadata"]["mass"] == pytest.approx(
        46.04, abs=0.01
    )


def test_parse_step_returns_no_step_for_leaf_node(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """A leaf node (molecule with no children) is a terminal starting material
    and should not generate a reaction step or any dependencies."""
    output = route_parser_factory().parse_step({"smiles": "CC"})

    assert output == {"dependencies": {}, "steps": []}


def test_parse_step_creates_product_only_step_for_empty_children(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """A node with an empty children list represents a reaction with no known
    precursors. It should still create a step with the product but leave
    reactants and reagents empty."""
    output = route_parser_factory().parse_step({"smiles": "CCO", "children": []})

    assert output["dependencies"] == {}
    assert len(output["steps"]) == 1
    assert output["steps"][0]["products"][0]["smiles"] == "CCO"
    assert output["steps"][0]["reactants"] == []
    assert output["steps"][0]["reagents"] == []


def test_parse_step_appends_to_existing_accumulators(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """When called with pre-existing step and dependency accumulators (as in
    recursive traversal), parse_step should mutate them in-place rather than
    creating new collections. Verifies the identity (`is`) contract and that
    the new step is correctly linked as a child dependency."""
    steps = [
        {
            "step": "1",
            "reactants": [],
            "reagents": [],
            "products": [{"smiles": "c1ccccc1", "product_metadata": {}}],
            "conditions": [],
            "reactionmetrics": [{"scalabilityindex": "", "closestliterature": ""}],
        }
    ]
    dependencies = {"1": []}

    output = route_parser_factory().parse_step(
        {"smiles": "CC(=O)O", "children": [{"children": [{"smiles": "CO"}]}]},
        step_list=steps,
        dependency_list=dependencies,
        parent_id=1,
    )

    assert output["steps"] is steps
    assert output["dependencies"] is dependencies
    assert output["dependencies"] == {"1": ["2"], "2": []}
    assert output["steps"][1]["step"] == "2"
    assert output["steps"][0]["reactants"][0]["smiles"] == "CC(=O)O"


def test_format_output_builds_dependency_chain_for_nested_route(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """A two-level route (propanol -> ethanol -> methanol) should produce two
    steps with step 1 depending on step 2, verifying the recursive traversal
    correctly chains multi-step retrosynthesis routes."""
    data = {
        "smiles": "CCCO",
        "children": [
            {
                "children": [
                    {
                        "smiles": "CCO",
                        "children": [{"children": [{"smiles": "CO"}]}],
                    }
                ]
            }
        ],
    }

    output = route_parser_factory().format_output(data)

    assert [step["products"][0]["smiles"] for step in output["steps"]] == [
        "CCCO",
        "CCO",
    ]
    assert output["dependencies"] == {"1": ["2"], "2": []}


def test_format_output_classifies_basic_molecules_as_reagents(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """Molecules in the basic_molecules set (e.g. water) should be placed in
    the reagents list, while other precursors (e.g. methanol) go into
    reactants. Tests the domain-specific reactant vs reagent classification."""
    data = {
        "smiles": "CC(=O)O",
        "children": [{"children": [{"smiles": "CO"}, {"smiles": "O"}]}],
    }

    output = route_parser_factory(basic_molecules={"O"}).format_output(data)
    step = output["steps"][0]

    assert [r["smiles"] for r in step["reactants"]] == ["CO"]
    assert [r["smiles"] for r in step["reagents"]] == ["O"]


def test_format_output_computes_scalability_only_from_reactants(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """The scalability calculator should only be invoked for reactants, not
    reagents (basic molecules like water). The stored scalabilityindex should
    reflect the last reactant processed, ignoring reagents entirely."""
    calls = []

    def record_scalability(reactant: str, product: str) -> int:
        calls.append((reactant, product))
        return 5

    data = {
        "smiles": "CC(=O)O",
        "children": [{"children": [{"smiles": "CO"}, {"smiles": "O"}]}],
    }

    output = route_parser_factory(
        scalability_calculator=record_scalability,
        basic_molecules={"O"},
    ).format_output(data)

    assert calls == [("CO", "CC(=O)O")]
    assert output["steps"][0]["reactionmetrics"][0]["scalabilityindex"] == 5


def test_fix_dependencies_ignores_unmatched_reactants() -> None:
    """Reactants whose SMILES don't match any step's product (e.g. butane
    not produced by any step) should be silently ignored in the dependency
    map rather than raising an error."""
    steps = [
        {"step": "1", "products": [{"smiles": "CO"}], "reactants": []},
        {
            "step": "2",
            "products": [{"smiles": "CCO"}],
            "reactants": [{"smiles": "CCCC"}],
        },
    ]

    assert fix_dependencies({}, steps) == {"1": [], "2": []}


def test_fix_dependencies_uses_last_duplicate_product_producer() -> None:
    """When multiple steps produce the same molecule (e.g. methanol in steps
    1 and 2), the dependency resolver links to the last producer. This
    documents the dict-overwrite semantics of product_to_step mapping."""
    steps = [
        {"step": "1", "products": [{"smiles": "CO"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "CO"}], "reactants": []},
        {"step": "3", "products": [{"smiles": "CCO"}], "reactants": [{"smiles": "CO"}]},
    ]

    assert fix_dependencies({}, steps)["3"] == ["2"]


def test_format_output_handles_reaction_wrapper_children(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """Route trees may contain intermediate reaction wrapper nodes (marked
    with is_reaction=True). The parser should traverse through these wrappers
    and correctly extract the precursor molecules underneath."""
    data = {
        "smiles": "CC(=O)O",
        "children": [{"is_reaction": True, "children": [{"smiles": "CO"}]}],
    }

    output = route_parser_factory().format_output(data)

    assert [r["smiles"] for r in output["steps"][0]["reactants"]] == ["CO"]


def test_format_output_preserves_duplicate_reactants(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """When the same molecule appears multiple times as a precursor (e.g. two
    equivalents of methanol), both should be preserved in the reactants list
    rather than being deduplicated."""
    data = {
        "smiles": "CC(=O)O",
        "children": [{"children": [{"smiles": "CO"}, {"smiles": "CO"}]}],
    }

    output = route_parser_factory().format_output(data)

    assert [r["smiles"] for r in output["steps"][0]["reactants"]] == ["CO", "CO"]


def test_format_output_raises_for_missing_smiles_anywhere_in_route(
    route_parser_factory: Callable[..., RetrosynthesisRouteParser],
) -> None:
    """A route node missing the required 'smiles' key should raise a
    ValueError, even if the root node is valid. This ensures malformed
    route trees from upstream are caught early."""
    with pytest.raises(ValueError, match="smiles"):
        route_parser_factory().format_output(
            {"smiles": "CC(=O)O", "children": [{"children": [{}]}]}
        )


def test_public_wrapper_functions_parse_routes_end_to_end() -> None:
    """The module-level parse_step() and format_output() compatibility
    wrappers should produce the same results as calling the class directly,
    ensuring the legacy API remains functional."""
    data = {
        "smiles": "CCO",
        "children": [{"children": [{"smiles": "CC"}]}],
    }

    parsed = parse_step(data)
    formatted = format_output(data)

    assert [step["products"][0]["smiles"] for step in parsed["steps"]] == ["CCO"]
    assert [step["products"][0]["smiles"] for step in formatted["steps"]] == ["CCO"]
    assert formatted["dependencies"] == {"1": []}
