"""Tests for package-level reaction metadata helpers."""

from __future__ import annotations

import pytest

from deepretro import metadata


def test_parse_reaction_smiles_supports_reaction_smiles_forms() -> None:
    participants = metadata.parse_reaction_smiles("CCO.CC(=O)O>>CCOC(C)=O")

    assert participants.reactants == [{"smiles": "CCO"}, {"smiles": "CC(=O)O"}]
    assert participants.reagents == []
    assert participants.product == [{"smiles": "CCOC(C)=O"}]

    participants_with_agent = metadata.parse_reaction_smiles("CCO>O>CC=O")

    assert participants_with_agent.reactants == [{"smiles": "CCO"}]
    assert participants_with_agent.reagents == [{"smiles": "O"}]
    assert participants_with_agent.product == [{"smiles": "CC=O"}]


def test_parse_reaction_smiles_rejects_invalid_shape() -> None:
    with pytest.raises(ValueError, match="reaction SMILES"):
        metadata.parse_reaction_smiles("CCO")


def test_build_reagent_records_filters_invalid_smiles_and_adds_metadata() -> None:
    records = metadata.build_reagent_records(["O", "not-smiles"])

    assert records == [
        {
            "smiles": "O",
            "reagent_metadata": {
                "name": "",
                "chemical_formula": "H2O",
                "mass": pytest.approx(18.01056468403),
            },
        }
    ]


def test_valid_conditions_payload_requires_all_condition_fields() -> None:
    assert metadata.valid_conditions_payload(
        {
            "temperature": "25 C",
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 h",
        }
    )
    assert not metadata.valid_conditions_payload({"temperature": "25 C"})
    assert not metadata.valid_conditions_payload("not a mapping")
    assert not metadata.valid_conditions_payload(
        {
            "temperature": "",
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 h",
        }
    )
    assert not metadata.valid_conditions_payload(
        {
            "temperature": None,
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 h",
        }
    )
    assert not metadata.valid_conditions_payload(
        {
            "temperature": [],
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 h",
        }
    )


def test_parse_literature_reaction_extracts_payload() -> None:
    response_text = (
        "{'literature_reaction': {'doi': '10.1000/example', 'year': 2024}, "
        "'explanation': 'example'}"
    )

    assert metadata.parse_literature_reaction(response_text) == {
        "doi": "10.1000/example",
        "year": 2024,
    }


def test_parse_metadata_response_accepts_json_fences_and_surrounding_text() -> None:
    assert metadata.parse_metadata_response('{"data": ["O"], "ok": true}') == {
        "data": ["O"],
        "ok": True,
    }
    assert metadata.parse_metadata_response(
        '```json\n{"data": ["O"], "ok": true}\n```'
    ) == {
        "data": ["O"],
        "ok": True,
    }
    assert metadata.parse_metadata_response(
        'Here is the payload:\n{"data": ["O"], "ok": true}\nDone.'
    ) == {
        "data": ["O"],
        "ok": True,
    }


def test_metadata_prompt_builders_include_reaction_context() -> None:
    reagent_messages = metadata.build_reagent_messages(["BrCBr"], "ClCCl")
    assert [message["role"] for message in reagent_messages] == ["system", "user"]
    assert "BrCBr" in reagent_messages[1]["content"]
    assert "ClCCl" in reagent_messages[1]["content"]
    assert "{reactants}" not in reagent_messages[1]["content"]
    assert "{product}" not in reagent_messages[1]["content"]

    condition_messages = metadata.build_conditions_messages(["BrCBr"], "ClCCl", ["N#N"])
    assert "BrCBr" in condition_messages[1]["content"]
    assert "ClCCl" in condition_messages[1]["content"]
    assert "N#N" in condition_messages[1]["content"]
    assert "{reactants}" not in condition_messages[1]["content"]
    assert "{product}" not in condition_messages[1]["content"]
    assert "{reagents}" not in condition_messages[1]["content"]

    literature_messages = metadata.build_literature_messages(
        ["BrCBr"],
        "ClCCl",
        ["N#N"],
        {
            "temperature": "37 C",
            "pressure": "1 atm",
            "solvent": "acetonitrile",
            "time": "2 h",
        },
    )
    assert "BrCBr" in literature_messages[1]["content"]
    assert "ClCCl" in literature_messages[1]["content"]
    assert "N#N" in literature_messages[1]["content"]
    assert "37 C" in literature_messages[1]["content"]
    assert "acetonitrile" in literature_messages[1]["content"]
    assert "{reactants}" not in literature_messages[1]["content"]
    assert "{product}" not in literature_messages[1]["content"]
    assert "{reagents}" not in literature_messages[1]["content"]
    assert "{conditions}" not in literature_messages[1]["content"]
