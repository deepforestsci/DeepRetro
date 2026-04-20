from deepretro.utils.parse import (
    RetrosynthesisRouteParser,
    parse_step,
    format_output,
    fix_dependencies,
)


def dummy_formula(smiles):
    return "X"

def dummy_mass(smiles):
    return 1.0

def dummy_scalability(a, b):
    return "N/A"


def test_simple_route_parsing():
    """Checks basic parsing of a simple route with two reactants."""
    data = {
        "smiles": "CCO",
        "children": [
            {
                "children": [
                    {"smiles": "CC"},
                    {"smiles": "O"},
                ]
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)

    assert "steps" in output
    assert "dependencies" in output

    assert len(output["steps"]) == 1
    step = output["steps"][0]

    assert step["products"][0]["smiles"] == "CCO"
    assert len(step["reactants"]) == 2

    product_metadata = step["products"][0]["product_metadata"]
    assert "chemical_formula" in product_metadata
    assert "mass" in product_metadata


def test_fix_dependencies_simple():
    """Checks dependency reconstruction from simple step relationships."""
    steps = [
        {"step": "1", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "B"}], "reactants": [{"smiles": "A"}]},
    ]

    deps = fix_dependencies({}, steps)

    assert deps["2"] == ["1"]


def test_leaf_node_no_step_created():
    """Ensures nodes without children do not produce steps."""
    data = {"smiles": "CC"}

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.parse_step(data)

    assert output["steps"] == []
    assert output["dependencies"] == {}


def test_empty_children_creates_step():
    """Checks that nodes with empty children still create a step."""
    data = {"smiles": "CCO", "children": []}

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.parse_step(data)

    assert len(output["steps"]) == 1

    step = output["steps"][0]
    assert step["products"][0]["smiles"] == "CCO"
    assert step["reactants"] == []
    assert step["reagents"] == []


def test_basic_molecule_goes_to_reagent():
    """Verifies basic molecules are classified as reagents."""
    data = {
        "smiles": "CCO",
        "children": [
            {
                "children": [
                    {"smiles": "H2O"},
                ]
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules={"H2O"},
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)
    step = output["steps"][0]

    assert len(step["reagents"]) == 1
    assert step["reagents"][0]["smiles"] == "H2O"


def test_missing_smiles_raises():
    """Ensures missing smiles field raises an error."""
    data = {"children": []}

    parser = RetrosynthesisRouteParser()

    try:
        parser.format_output(data)
    except ValueError as e:
        assert "smiles" in str(e)
    else:
        assert False


def test_multi_step_route():
    """Checks recursive parsing generates multiple steps."""
    data = {
        "smiles": "C",
        "children": [
            {
                "children": [
                    {
                        "smiles": "B",
                        "children": [
                            {
                                "children": [
                                    {"smiles": "A"}
                                ]
                            }
                        ],
                    }
                ]
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)

    assert len(output["steps"]) == 2


def test_parse_step_wrapper():
    """Checks wrapper parse_step returns structured output."""
    data = {
        "smiles": "CCO",
        "children": [{"children": [{"smiles": "CC"}]}],
    }

    result = parse_step(data)

    assert "steps" in result


def test_format_output_wrapper():
    """Checks wrapper format_output returns dependencies."""
    data = {
        "smiles": "CCO",
        "children": [{"children": [{"smiles": "CC"}]}],
    }

    result = format_output(data)

    assert "dependencies" in result


def test_one_parent_three_children():
    """Checks parsing of a single step with three reactants."""
    data = {
        "smiles": "P",
        "children": [
            {
                "children": [
                    {"smiles": "A"},
                    {"smiles": "B"},
                    {"smiles": "C"},
                ]
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)
    step = output["steps"][0]

    assert len(step["reactants"]) == 3
    assert {r["smiles"] for r in step["reactants"]} == {"A", "B", "C"}


def test_mixed_reactant_and_reagent():
    """Checks correct split between reactants and reagents."""
    data = {
        "smiles": "P",
        "children": [
            {
                "children": [
                    {"smiles": "A"},
                    {"smiles": "H2O"},
                ]
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules={"H2O"},
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)
    step = output["steps"][0]

    assert len(step["reactants"]) == 1
    assert len(step["reagents"]) == 1


def test_dependency_chain_three_steps():
    """Checks dependency reconstruction across multiple linked steps."""
    steps = [
        {"step": "1", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "B"}], "reactants": [{"smiles": "A"}]},
        {"step": "3", "products": [{"smiles": "C"}], "reactants": [{"smiles": "B"}]},
    ]

    deps = fix_dependencies({}, steps)

    assert deps["2"] == ["1"]
    assert deps["3"] == ["2"]


def test_duplicate_reactants_handled():
    """Checks behavior when same reactant appears multiple times."""
    data = {
        "smiles": "P",
        "children": [
            {
                "children": [
                    {"smiles": "A"},
                    {"smiles": "A"},
                ]
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)
    step = output["steps"][0]

    assert len(step["reactants"]) == 2


def test_no_matching_dependency():
    """Checks that unmatched reactants do not create dependencies."""
    steps = [
        {"step": "1", "products": [{"smiles": "A"}], "reactants": []},
        {"step": "2", "products": [{"smiles": "B"}], "reactants": [{"smiles": "X"}]},
    ]

    deps = fix_dependencies({}, steps)

    assert deps["2"] == []


def test_single_node_with_reaction_flag():
    """Checks that reaction nodes do not attach as reactants."""
    data = {
        "smiles": "P",
        "children": [
            {
                "is_reaction": True,
                "children": [{"smiles": "A"}],
            }
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)

    step = output["steps"][0]
    assert len(step["reactants"]) == 1


def test_multiple_independent_steps():
    """Checks parsing when multiple independent branches exist."""
    data = {
        "smiles": "P",
        "children": [
            {"children": [{"smiles": "A"}]},
            {"children": [{"smiles": "B"}]},
        ],
    }

    parser = RetrosynthesisRouteParser(
        basic_molecules=set(),
        chemical_formula_calculator=dummy_formula,
        mass_calculator=dummy_mass,
        scalability_calculator=dummy_scalability,
    )

    output = parser.format_output(data)

    assert len(output["steps"]) >= 1


def test_missing_children_key_in_nested_node():
    """Checks behavior when nested node lacks children key."""
    data = {
        "smiles": "P",
        "children": [
            {
                "children": [{}]
            }
        ],
    }

    parser = RetrosynthesisRouteParser()

    try:
        parser.format_output(data)
    except ValueError:
        assert True
    else:
        assert False