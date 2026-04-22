import pytest

from deepretro.utils.parse import (
    RetrosynthesisRouteParser,
    fix_dependencies,
    format_output,
    parse_step,
)


def build_parser(scalability_calculator=None, basic_molecules=None):
    return RetrosynthesisRouteParser(
        basic_molecules=set() if basic_molecules is None else basic_molecules,
        chemical_formula_calculator=lambda smiles: "X",
        mass_calculator=lambda smiles: 1.0,
        scalability_calculator=(
            (lambda reactant, product: "N/A")
            if scalability_calculator is None
            else scalability_calculator
        ),
    )


def test_format_output_parses_single_step_with_metadata():
    data = {
        "smiles": "CCO",
        "children": [{"children": [{"smiles": "CC"}, {"smiles": "O"}]}],
    }

    output = build_parser().format_output(data)

    assert output["dependencies"] == {"1": []}
    assert len(output["steps"]) == 1

    step = output["steps"][0]
    assert step["step"] == "1"
    assert [product["smiles"] for product in step["products"]] == ["CCO"]
    assert [reactant["smiles"] for reactant in step["reactants"]] == ["CC", "O"]
    assert step["reagents"] == []
    assert step["reactionmetrics"][0]["scalabilityindex"] == "N/A"
    assert step["products"][0]["product_metadata"]["chemical_formula"] == "X"
    assert step["products"][0]["product_metadata"]["mass"] == 1.0


def test_parse_step_returns_no_step_for_leaf_node():
    output = build_parser().parse_step({"smiles": "CC"})

    assert output == {"dependencies": {}, "steps": []}


def test_parse_step_creates_product_only_step_for_empty_children():
    output = build_parser().parse_step({"smiles": "CCO", "children": []})

    assert output["dependencies"] == {}
    assert len(output["steps"]) == 1
    assert output["steps"][0]["products"][0]["smiles"] == "CCO"
    assert output["steps"][0]["reactants"] == []
    assert output["steps"][0]["reagents"] == []


def test_parse_step_appends_to_existing_accumulators():
    steps = [
        {
            "step": "1",
            "reactants": [],
            "reagents": [],
            "products": [{"smiles": "seed", "product_metadata": {}}],
            "conditions": [],
            "reactionmetrics": [{"scalabilityindex": "", "closestliterature": ""}],
        }
    ]
    dependencies = {"1": []}

    output = build_parser().parse_step(
        {"smiles": "P", "children": [{"children": [{"smiles": "A"}]}]},
        step_list=steps,
        dependency_list=dependencies,
        parent_id=1,
    )

    assert output["steps"] is steps
    assert output["dependencies"] is dependencies
    assert output["dependencies"] == {"1": ["2"], "2": []}
    assert output["steps"][1]["step"] == "2"
    assert output["steps"][0]["reactants"][0]["smiles"] == "P"


def test_format_output_builds_dependency_chain_for_nested_route():
    data = {
        "smiles": "C",
        "children": [
            {
                "children": [
                    {
                        "smiles": "B",
                        "children": [{"children": [{"smiles": "A"}]}],
                    }
                ]
            }
        ],
    }

    output = build_parser().format_output(data)

    assert [step["products"][0]["smiles"] for step in output["steps"]] == ["C", "B"]
    assert output["dependencies"] == {"1": ["2"], "2": []}


def test_format_output_classifies_basic_molecules_as_reagents():
    data = {
        "smiles": "P",
        "children": [{"children": [{"smiles": "A"}, {"smiles": "H2O"}]}],
    }

    output = build_parser(basic_molecules={"H2O"}).format_output(data)
    step = output["steps"][0]

    assert [reactant["smiles"] for reactant in step["reactants"]] == ["A"]
    assert [reagent["smiles"] for reagent in step["reagents"]] == ["H2O"]


def test_format_output_uses_reactants_not_reagents_for_scalability():
    calls = []

    def record_scalability(reactant, product):
        calls.append((reactant, product))
        return f"{reactant}->{product}"

    data = {
        "smiles": "P",
        "children": [{"children": [{"smiles": "A"}, {"smiles": "H2O"}]}],
    }

    output = build_parser(
        scalability_calculator=record_scalability,
        basic_molecules={"H2O"},
    ).format_output(data)

    assert calls == [("A", "P")]
    assert output["steps"][0]["reactionmetrics"][0]["scalabilityindex"] == "A->P"


def test_fix_dependencies_ignores_unmatched_reactants():
    steps = [
        {"step": "1", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "B"}], "reactants": [{"smiles": "X"}]},
    ]

    assert fix_dependencies({}, steps) == {"1": [], "2": []}


def test_fix_dependencies_preserves_all_duplicate_producers():
    steps = [
        {"step": "1", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "3", "products": [{"smiles": "B"}], "reactants": [{"smiles": "A"}]},
    ]

    assert fix_dependencies({}, steps)["3"] == ["1", "2"]


def test_format_output_handles_reaction_wrapper_children():
    data = {
        "smiles": "P",
        "children": [{"is_reaction": True, "children": [{"smiles": "A"}]}],
    }

    output = build_parser().format_output(data)

    assert [reactant["smiles"] for reactant in output["steps"][0]["reactants"]] == ["A"]


def test_format_output_preserves_duplicate_reactants():
    data = {
        "smiles": "P",
        "children": [{"children": [{"smiles": "A"}, {"smiles": "A"}]}],
    }

    output = build_parser().format_output(data)

    assert [reactant["smiles"] for reactant in output["steps"][0]["reactants"]] == [
        "A",
        "A",
    ]


def test_format_output_raises_for_missing_smiles_anywhere_in_route():
    with pytest.raises(ValueError, match="smiles"):
        build_parser().format_output(
            {"smiles": "P", "children": [{"children": [{}]}]}
        )


def test_public_wrapper_functions_parse_routes_end_to_end():
    data = {
        "smiles": "CCO",
        "children": [{"children": [{"smiles": "CC"}]}],
    }

    parsed = parse_step(data)
    formatted = format_output(data)

    assert [step["products"][0]["smiles"] for step in parsed["steps"]] == ["CCO"]
    assert [step["products"][0]["smiles"] for step in formatted["steps"]] == ["CCO"]
    assert formatted["dependencies"] == {"1": []}
