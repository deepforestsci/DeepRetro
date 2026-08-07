"""Tests for :mod:`deepretro.score`."""

import pytest

from deepretro import score


def test_canonicalize_smiles_normalizes_valid_input():
    assert score.canonicalize_smiles("C1=CC=CC=C1") == "c1ccccc1"
    assert score.canonicalize_smiles(" CCO ") == "CCO"


@pytest.mark.parametrize("smiles", ["", "   ", "not_a_smiles!!!", None])
def test_canonicalize_smiles_rejects_invalid_input(smiles):
    assert score.canonicalize_smiles(smiles) is None


def test_sa_score_uses_the_real_rdkit_backend():
    value = score.sa_score("CCO")

    assert isinstance(value, float)
    assert 1.0 <= value <= 10.0
    assert score.sa_score("not_a_smiles!!!") is None


def test_sc_score_handles_an_unconfigured_optional_backend():
    value = score.sc_score("CCO")

    if value is None:
        molecule = score.score_molecule("CCO")
        assert any(
            warning.startswith("SCScore unavailable")
            for warning in molecule["warnings"]
        )
    else:
        assert 1.0 <= value <= 5.0

    assert score.sc_score("not_a_smiles!!!") is None


def test_score_molecule_end_to_end():
    result = score.score_molecule(" CCO ")

    assert result["input_smiles"] == " CCO "
    assert result["canonical_smiles"] == "CCO"
    assert result["valid"] is True
    assert isinstance(result["sa_score"], float)
    assert result["sc_score"] is None or isinstance(result["sc_score"], float)
    assert result["errors"] == []


def test_score_molecule_reports_invalid_smiles():
    result = score.score_molecule("not_a_smiles!!!")

    assert result["canonical_smiles"] is None
    assert result["valid"] is False
    assert result["sa_score"] is None
    assert result["sc_score"] is None
    assert result["errors"] == ["Invalid SMILES"]


def test_score_step_aggregates_real_molecule_scores():
    result = score.score_step("CCO", "CC=O.O")

    product_sa = result["product"]["sa_score"]
    reactant_sa = [reactant["sa_score"] for reactant in result["reactants"]]

    assert result["errors"] == []
    assert result["product"]["canonical_smiles"] == "CCO"
    assert [item["canonical_smiles"] for item in result["reactants"]] == [
        "CC=O",
        "O",
    ]
    assert len(result["reactants"]) == 2
    assert all(isinstance(value, float) for value in [product_sa, *reactant_sa])
    assert result["mean_reactant_sa"] == pytest.approx(sum(reactant_sa) / 2)
    assert result["max_reactant_sa"] == max(reactant_sa)
    assert result["min_reactant_sa"] == min(reactant_sa)
    assert result["delta_sa_mean"] == pytest.approx(
        product_sa - result["mean_reactant_sa"]
    )
    assert result["delta_sa_max"] == pytest.approx(
        product_sa - result["max_reactant_sa"]
    )
    assert result["sa_simplifies"] is (result["delta_sa_mean"] > 0)

    reactant_sc = [reactant["sc_score"] for reactant in result["reactants"]]
    if result["product"]["sc_score"] is None or any(
        value is None for value in reactant_sc
    ):
        assert result["mean_reactant_sc"] is None
        assert result["delta_sc_mean"] is None
        assert result["sc_simplifies"] is False


def test_score_step_accepts_list_and_dot_separated_reactants():
    from_string = score.score_step("CCO", "CC=O.O")
    from_list = score.score_step("CCO", ["CC=O", "O"])

    assert from_list["product"] == from_string["product"]
    assert from_list["reactants"] == from_string["reactants"]
    assert from_list["mean_reactant_sa"] == from_string["mean_reactant_sa"]
    assert from_list["delta_sa_mean"] == from_string["delta_sa_mean"]


def test_score_step_reports_missing_reactants():
    result = score.score_step("CCO", "")

    assert result["product"] is None
    assert result["reactants"] == []
    assert result["errors"] == ["Malformed step: missing reactant SMILES"]


def test_score_step_keeps_valid_reactants_when_product_is_invalid():
    result = score.score_step("not_a_smiles!!!", "CCO")

    assert result["product"]["valid"] is False
    assert result["reactants"][0]["valid"] is True
    assert result["delta_sa_mean"] is None
    assert result["sa_simplifies"] is False
    assert result["errors"] == ["Invalid SMILES"]


def test_score_pathway_accepts_simple_steps_end_to_end():
    pathway = [
        {"product": "CCO", "reactants": "CC=O.O"},
        {"product_smiles": "CC=O", "reactants_smiles": ["CC#N", "O"]},
    ]

    result = score.score_pathway(pathway)

    assert result["errors"] == []
    assert len(result["step_scores"]) == 2
    assert result["route_summary"]["n_steps"] == 2
    assert all(step["product"]["valid"] for step in result["step_scores"])


def test_score_pathway_accepts_viewer_steps_end_to_end():
    pathway = {
        "steps": [
            {
                "products": [{"smiles": "CCO"}],
                "reactants": [{"smiles": "CC=O"}, {"smiles": "O"}],
            }
        ]
    }

    result = score.score_pathway(pathway)

    step = result["step_scores"][0]
    assert result["errors"] == []
    assert result["route_summary"]["n_steps"] == 1
    assert step["product"]["canonical_smiles"] == "CCO"
    assert [item["canonical_smiles"] for item in step["reactants"]] == [
        "CC=O",
        "O",
    ]


@pytest.mark.parametrize(
    ("step", "error"),
    [
        (
            {"products": [], "reactants": [{"smiles": "CC=O"}]},
            "Malformed step 0: missing product SMILES",
        ),
        (
            {"products": [{"smiles": "CCO"}], "reactants": []},
            "Malformed step 0: missing reactant SMILES",
        ),
        (
            {
                "products": [{"smiles": "CCO"}],
                "reactants": [{"smiles": "CC=O"}],
                "product": "CCO",
            },
            "Malformed step 0: mixed viewer-style and simple-style step",
        ),
    ],
)
def test_score_pathway_reports_invalid_viewer_steps(step, error):
    result = score.score_pathway({"steps": [step]})

    assert result["step_scores"][0]["product"] is None
    assert result["errors"] == [error]


@pytest.mark.parametrize(
    ("pathway", "error"),
    [
        ("not a pathway", "Malformed pathway: expected dict or list"),
        ({}, "Malformed pathway: missing steps"),
        ({"steps": {}}, "Malformed pathway: steps must be a list"),
    ],
)
def test_score_pathway_rejects_malformed_input(pathway, error):
    result = score.score_pathway(pathway)

    assert result["step_scores"] == []
    assert result["route_summary"]["n_steps"] == 0
    assert result["errors"] == [error]


def test_route_summary_aggregates_step_scores():
    step_scores = [
        {
            "product": {"sa_score": 6.0, "sc_score": 4.0},
            "delta_sa_mean": 3.0,
            "delta_sc_mean": 2.0,
            "sa_simplifies": True,
            "sc_simplifies": True,
        },
        {
            "product": {"sa_score": 2.0, "sc_score": 1.0},
            "delta_sa_mean": -2.0,
            "delta_sc_mean": -2.0,
            "sa_simplifies": False,
            "sc_simplifies": False,
        },
    ]

    assert score._summarize_steps(step_scores) == {
        "n_steps": 2,
        "mean_product_sa": 4.0,
        "mean_product_sc": 2.5,
        "mean_delta_sa": 0.5,
        "mean_delta_sc": 0.0,
        "fraction_steps_sa_simplifying": 0.5,
        "fraction_steps_sc_simplifying": 0.5,
        "complexity_increasing_steps": [1],
    }


def test_empty_pathway_scores_uses_the_standard_schema():
    result = score.empty_pathway_scores("scoring failed")

    assert set(result) == {"step_scores", "route_summary", "warnings", "errors"}
    assert result["step_scores"] == []
    assert result["route_summary"]["n_steps"] == 0
    assert result["warnings"] == []
    assert result["errors"] == ["scoring failed"]
