"""Tests for package-level reaction metadata helpers."""

from __future__ import annotations

import os
from typing import cast

import pytest

from deepretro import metadata
from deepretro.utils.cache import CacheManager, make_cache_key

OPUS_MODEL = os.getenv("DEEPRETRO_METADATA_TEST_MODEL", metadata.DEFAULT_METADATA_MODEL)


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


def test_recommend_reaction_metadata_from_reaction_string_uses_all_three_agents(
    spy_reagent_recommender: metadata.ReagentRecommender,
    spy_conditions_recommender: metadata.ConditionsRecommender,
    spy_literature_recommender: metadata.LiteratureRecommender,
    recommender_calls: list[str],
) -> None:
    status, recommendation = metadata.recommend_reaction_metadata(
        "CCO.CC(=O)O>>CCOC(C)=O",
        model="test-model",
        temperature=0.2,
        reagent_recommender=spy_reagent_recommender,
        conditions_recommender=spy_conditions_recommender,
        literature_recommender=spy_literature_recommender,
        cache=None,
    )

    assert status == 200
    assert recommendation == {
        "reaction_smiles": "CCO.CC(=O)O>>CCOC(C)=O",
        "reactants": [{"smiles": "CCO"}, {"smiles": "CC(=O)O"}],
        "product": [{"smiles": "CCOC(C)=O"}],
        "reagents": [{"smiles": "O", "reagent_metadata": {"name": ""}}],
        "conditions": {
            "temperature": "25 C",
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 h",
        },
        "literature": {"doi": "10.1000/example"},
    }
    assert recommender_calls == [
        "reagent:test-model:0.2",
        "conditions:test-model:0.2",
        "literature:test-model:0.2",
    ]


def test_recommend_reaction_metadata_returns_failure_status(
    failing_reagent_recommender: metadata.ReagentRecommender,
) -> None:
    status, recommendation = metadata.recommend_reaction_metadata(
        "CCO>>CC=O",
        reagent_recommender=failing_reagent_recommender,
        cache=None,
    )

    assert status == 404
    assert recommendation == {
        "stage": "reagents",
        "error": "metadata recommendation failed",
    }


def test_recommend_reaction_metadata_returns_seeded_cache_hit() -> None:
    cache = CacheManager()
    reaction_smiles = "CCO>>CC=O"
    expected: metadata.MetadataRecommendation = {
        "reaction_smiles": reaction_smiles,
        "reactants": [{"smiles": "cached-reactant"}],
        "product": [{"smiles": "cached-product"}],
        "reagents": [{"smiles": "cached-reagent"}],
        "conditions": {
            "temperature": "cached-temperature",
            "pressure": "cached-pressure",
            "solvent": "cached-solvent",
            "time": "cached-time",
        },
        "literature": {"doi": "cached-doi"},
    }
    cache_key = make_cache_key(
        "recommend_reaction_metadata",
        reaction_smiles,
        model=metadata.DEFAULT_METADATA_MODEL,
        temperature=0.0,
    )
    cache.set(cache_key, (200, expected), tag="metadata")

    status, recommendation = metadata.recommend_reaction_metadata(
        reaction_smiles,
        cache=cache,
    )

    assert status == 200
    assert recommendation == expected


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


def test_metadata_agents_hit_real_opus() -> None:
    if not os.getenv("ANTHROPIC_API_KEY"):
        pytest.skip("ANTHROPIC_API_KEY not set — skipping live Opus metadata test")

    reactants: list[metadata.MoleculeRecord] = [
        {"smiles": "CCO"},
        {"smiles": "CC(=O)O"},
    ]
    product: list[metadata.MoleculeRecord] = [{"smiles": "CCOC(C)=O"}]

    reagent_status, reagents = metadata.reagent_agent(
        reactants=reactants,
        product=product,
        model=OPUS_MODEL,
        temperature=0.0,
    )

    assert reagent_status == 200
    assert isinstance(reagents, list)
    assert reagents
    assert all("smiles" in reagent for reagent in reagents)
    assert all("reagent_metadata" in reagent for reagent in reagents)

    status, conditions = metadata.conditions_agent(
        reactants=reactants,
        product=product,
        reagents=reagents,
        model=OPUS_MODEL,
        temperature=0.0,
    )

    assert status == 200
    assert metadata.valid_conditions_payload(conditions)
    assert str(conditions["temperature"]).strip()
    assert str(conditions["pressure"]).strip()
    assert str(conditions["solvent"]).strip()
    assert str(conditions["time"]).strip()

    literature_status, literature = metadata.literature_agent(
        reactants=reactants,
        product=product,
        reagents=reagents,
        conditions=conditions,
        model=OPUS_MODEL,
        temperature=0.0,
    )

    assert literature_status == 200
    assert str(literature).strip()

    recommendation_status, recommendation = metadata.recommend_reaction_metadata(
        "CCO.CC(=O)O>>CCOC(C)=O",
        model=OPUS_MODEL,
        temperature=0.0,
        cache=None,
    )

    assert recommendation_status == 200
    recommendation_payload = cast(metadata.MetadataRecommendation, recommendation)
    assert metadata.valid_conditions_payload(recommendation_payload["conditions"])
    assert recommendation_payload["reagents"]
    assert str(recommendation_payload["literature"]).strip()
