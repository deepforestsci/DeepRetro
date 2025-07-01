import pytest
from unittest import mock

# Assuming src.utils.parse is the module to test
from src.utils import parse

# Test functions will be added here 

def test_fix_dependencies_linear_pathway():
    """Test fix_dependencies with a simple linear pathway A -> B -> C."""
    # Step C: Product C, Reactant B
    # Step B: Product B, Reactant A
    # Step A: Product A, (no reactants needed for this test structure)
    step_list = [
        {
            "step": "1", # Produces A
            "products": [{"smiles": "A"}],
            "reactants": [] 
        },
        {
            "step": "2", # Produces B from A
            "products": [{"smiles": "B"}],
            "reactants": [{"smiles": "A"}]
        },
        {
            "step": "3", # Produces C from B
            "products": [{"smiles": "C"}],
            "reactants": [{"smiles": "B"}]
        }
    ]
    initial_dependencies = {} # The input `dependencies` dict can be empty

    expected_fixed_dependencies = {
        "1": [],
        "2": ["1"],  # Step 2 depends on Step 1 (A)
        "3": ["2"]   # Step 3 depends on Step 2 (B)
    }

    result = parse.fix_dependencies(initial_dependencies, step_list)
    assert result == expected_fixed_dependencies

def test_fix_dependencies_branched_pathway():
    """Test fix_dependencies with a branched pathway (A & B -> C)."""
    # Step C: Product C, Reactants A, B
    # Step A: Product A
    # Step B: Product B
    step_list = [
        {
            "step": "1", # Produces A
            "products": [{"smiles": "A"}],
            "reactants": []
        },
        {
            "step": "2", # Produces B
            "products": [{"smiles": "B"}],
            "reactants": []
        },
        {
            "step": "3", # Produces C from A and B
            "products": [{"smiles": "C"}],
            "reactants": [{"smiles": "A"}, {"smiles": "B"}]
        }
    ]
    initial_dependencies = {}

    expected_fixed_dependencies = {
        "1": [],
        "2": [],
        "3": ["1", "2"] # Step 3 depends on Step 1 (A) and Step 2 (B)
    }
    # The order of dependencies for step "3" might vary, so check as a set
    result = parse.fix_dependencies(initial_dependencies, step_list)
    assert result["1"] == expected_fixed_dependencies["1"]
    assert result["2"] == expected_fixed_dependencies["2"]
    assert sorted(result["3"]) == sorted(expected_fixed_dependencies["3"])
    assert len(result) == len(expected_fixed_dependencies)

def test_fix_dependencies_no_dependencies():
    """Test fix_dependencies when steps have no inter-dependencies (e.g., only starting materials)."""
    step_list = [
        {
            "step": "1",
            "products": [{"smiles": "A"}],
            "reactants": [{"smiles": "X"}] # X is not a product of another step
        },
        {
            "step": "2",
            "products": [{"smiles": "B"}],
            "reactants": [{"smiles": "Y"}] # Y is not a product of another step
        }
    ]
    initial_dependencies = {}

    expected_fixed_dependencies = {
        "1": [],
        "2": []
    }

    result = parse.fix_dependencies(initial_dependencies, step_list)
    assert result == expected_fixed_dependencies

def test_fix_dependencies_empty_step_list():
    """Test fix_dependencies with an empty step_list."""
    step_list = []
    initial_dependencies = {}
    expected_fixed_dependencies = {}
    result = parse.fix_dependencies(initial_dependencies, step_list)
    assert result == expected_fixed_dependencies

def test_fix_dependencies_input_dependencies_ignored():
    """Test that the input `dependencies` argument is effectively ignored and overwritten."""
    step_list = [
        {
            "step": "1", 
            "products": [{"smiles": "A"}],
            "reactants": [] 
        },
        {
            "step": "2", 
            "products": [{"smiles": "B"}],
            "reactants": [{"smiles": "A"}]
        }
    ]
    # Provide some initial dependencies that should be cleared/rebuilt
    initial_dependencies = {"1": ["X"], "2": ["Y"], "random": ["Z"]}

    expected_fixed_dependencies = {
        "1": [],
        "2": ["1"]
    }
    result = parse.fix_dependencies(initial_dependencies, step_list)
    assert result == expected_fixed_dependencies 