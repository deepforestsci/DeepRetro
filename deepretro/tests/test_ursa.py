from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
from ursa import RetrosyntheticPath
from ursa.collapsing.path_collapser import PathCollapser
from ursa.validation.path_consistency_checker import PathConsistencyChecker

from deepretro import score
from deepretro.ursa import (
    UrsaAdapter,
    UrsaConfig,
    UrsaConfigurationError,
    UrsaRouteError,
    _to_retrocast_route,
)


LINEAR_PATHWAY = [
    {"product": "CCO", "reactants": ["CC=O", "O"]},
    {"product": "CC=O", "reactants": ["CC"]},
]

RELEASED_EXAMPLE_PATHWAY = [
    {
        "product": "CN1[C@H](CNc2ccccc2Nc2ccc(-n3nccn3)nc2)CCS1(=O)=O",
        "reactants": [
            "Nc1ccccc1Nc1ccc(-n2nccn2)nc1",
            "CN1[C@H](C=O)CCS1(=O)=O",
        ],
    },
    {
        "product": "Nc1ccccc1Nc1ccc(-n2nccn2)nc1",
        "reactants": ["O=[N+]([O-])c1ccccc1Nc1ccc(-n2nccn2)nc1"],
    },
    {
        "product": "O=[N+]([O-])c1ccccc1Nc1ccc(-n2nccn2)nc1",
        "reactants": ["Nc1ccc(-n2nccn2)nc1", "O=[N+]([O-])c1ccccc1F"],
    },
    {
        "product": "CN1[C@H](C=O)CCS1(=O)=O",
        "reactants": ["CN1[C@H](CO)CCS1(=O)=O"],
    },
    {
        "product": "CN1[C@H](CO)CCS1(=O)=O",
        "reactants": ["CN1[C@H](COCc2ccccc2)CCS1(=O)=O"],
    },
    {
        "product": "CN1[C@H](COCc2ccccc2)CCS1(=O)=O",
        "reactants": ["CI", "O=S1(=O)CC[C@@H](COCc2ccccc2)N1"],
    },
    {
        "product": "O=S1(=O)CC[C@@H](COCc2ccccc2)N1",
        "reactants": ["CS(=O)(=O)Cl", "N[C@@H](CO)COCc1ccccc1"],
    },
]

RELEASED_EXAMPLE_ROOT_REACTION = (
    "Nc1ccccc1Nc1ccc(-n2nccn2)nc1.CN1[C@H](C=O)CCS1(=O)=O"
    ">>CN1[C@H](CNc2ccccc2Nc2ccc(-n3nccn3)nc2)CCS1(=O)=O"
)


def existing_asset_config(tmp_path: Path) -> UrsaConfig:
    catalog = tmp_path / "building-blocks.csv"
    database = tmp_path / "chemcensor.sqlite"
    catalog.touch()
    database.touch()
    return UrsaConfig(
        enabled=True,
        bb_catalog_path=catalog,
        chemcensor_db_path=database,
    )


def configured_real_assets() -> UrsaConfig:
    catalog = os.environ.get("URSA_BB_CATALOG_PATH")
    database = os.environ.get("URSA_CHEMCENSOR_DB_PATH")
    if not catalog or not database:
        pytest.skip(
            "set URSA_BB_CATALOG_PATH and URSA_CHEMCENSOR_DB_PATH to run "
            "the real URSA integration test"
        )
    return UrsaConfig(
        enabled=True,
        bb_catalog_path=catalog,
        chemcensor_db_path=database,
    )


def test_score_import_eagerly_imports_official_ursa() -> None:
    probe = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import deepretro.score; "
            "assert 'deepretro.ursa' in sys.modules; assert 'ursa' in sys.modules",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert probe.returncode == 0, probe.stderr


def test_disabled_config_needs_no_asset_paths() -> None:
    assert UrsaConfig() == UrsaConfig(
        enabled=False,
        bb_catalog_path=None,
        chemcensor_db_path=None,
    )


@pytest.mark.parametrize(
    "pathway",
    [
        [{"product": "CCO", "reactants": ["CC=O", "O"]}],
        [],
        [{"product": "CCO"}],
        [{"product": "not-smiles", "reactants": ["CC"]}],
        [{"product": "CCO", "reactants": ["not-smiles"]}],
    ],
)
def test_score_pathway_is_unchanged_when_ursa_is_disabled(pathway) -> None:
    baseline = score.score_pathway(pathway)

    result = score.score_pathway(pathway, ursa_config=UrsaConfig())

    assert result == baseline
    assert "ursa" not in result


@pytest.mark.parametrize(
    "pathway",
    [
        [],
        [{"product": "CCO"}],
        [{"product": "not-smiles", "reactants": ["CC"]}],
        [{"product": "CCO", "reactants": ["not-smiles"]}],
    ],
)
def test_enabled_ursa_is_skipped_for_unscorable_pathways(
    tmp_path: Path,
    pathway,
) -> None:
    result = score.score_pathway(
        pathway,
        ursa_config=existing_asset_config(tmp_path),
    )

    assert "ursa" not in result


def test_disabled_config_does_not_validate_supplied_paths(tmp_path: Path) -> None:
    config = UrsaConfig(
        bb_catalog_path=tmp_path / "missing.csv",
        chemcensor_db_path=tmp_path / "missing.sqlite",
    )

    assert config.bb_catalog_path == (tmp_path / "missing.csv").resolve()
    assert config.chemcensor_db_path == (tmp_path / "missing.sqlite").resolve()


def test_enabled_config_normalizes_existing_asset_paths(tmp_path: Path) -> None:
    config = existing_asset_config(tmp_path)

    assert config.bb_catalog_path == (tmp_path / "building-blocks.csv").resolve()
    assert config.chemcensor_db_path == (tmp_path / "chemcensor.sqlite").resolve()


@pytest.mark.parametrize(
    ("catalog_exists", "database_exists", "message"),
    [
        (False, True, "building-block catalog"),
        (True, False, "ChemCensor database"),
    ],
)
def test_enabled_config_rejects_missing_assets(
    tmp_path: Path,
    catalog_exists: bool,
    database_exists: bool,
    message: str,
) -> None:
    catalog = tmp_path / "building-blocks.csv"
    database = tmp_path / "chemcensor.sqlite"
    if catalog_exists:
        catalog.touch()
    if database_exists:
        database.touch()

    with pytest.raises(UrsaConfigurationError, match=message):
        UrsaConfig(
            enabled=True,
            bb_catalog_path=catalog,
            chemcensor_db_path=database,
        )


def test_disabled_adapter_rejects_annotation() -> None:
    with pytest.raises(UrsaConfigurationError, match="disabled"):
        UrsaAdapter(UrsaConfig()).annotate(LINEAR_PATHWAY)


def test_single_viewer_step_uses_retrocast_shape() -> None:
    route = _to_retrocast_route(
        {
            "steps": [
                {
                    "products": [{"smiles": "CCO"}],
                    "reactants": [{"smiles": "CC"}, {"smiles": "O"}],
                }
            ]
        }
    )

    assert route == {
        "schema_version": "2",
        "target": {
            "smiles": "CCO",
            "product_of": {
                "reactants": [
                    {"smiles": "CC", "product_of": None},
                    {"smiles": "O", "product_of": None},
                ]
            },
        },
    }


def test_out_of_order_steps_link_canonically_equivalent_intermediate() -> None:
    route = _to_retrocast_route(
        [
            {"product": "CC=O", "reactants": ["CC"]},
            {"product": "CCO", "reactants": ["C(C)=O", "O"]},
        ]
    )

    intermediate = route["target"]["product_of"]["reactants"][0]
    assert intermediate == {
        "smiles": "CC=O",
        "product_of": {
            "reactants": [{"smiles": "CC", "product_of": None}],
        },
    }


def test_conversion_builds_a_branching_tree() -> None:
    official_path = RetrosyntheticPath.from_dict(
        _to_retrocast_route(
            [
                {"product": "CCOC", "reactants": ["CC=O", "CO"]},
                {"product": "CO", "reactants": ["C", "O"]},
                {"product": "CC=O", "reactants": ["CC"]},
            ]
        ),
        "branching-route",
    )

    assert official_path.num_steps == 3
    assert {node.smiles for node in official_path.root.children} == {"CC=O", "CO"}
    assert {node.smiles for node in official_path.get_starting_materials()} == {
        "C",
        "CC",
        "O",
    }


def test_conversion_preserves_stereochemistry_and_components() -> None:
    route = _to_retrocast_route(
        [
            {
                "product": "C[C@H](O)Cl",
                "reactants": ["C[C@H](O)Br", "[Na+]"],
            }
        ]
    )

    target = route["target"]
    assert target["smiles"] == "C[C@H](O)Cl"
    assert [node["smiles"] for node in target["product_of"]["reactants"]] == [
        "C[C@H](O)Br",
        "[Na+]",
    ]


def test_conversion_builds_a_real_official_path() -> None:
    official_path = RetrosyntheticPath.from_dict(
        _to_retrocast_route(LINEAR_PATHWAY),
        "real-conversion-smoke",
    )

    assert isinstance(official_path, RetrosyntheticPath)
    assert official_path.path_id == "real-conversion-smoke"
    assert official_path.root.smiles == "CCO"
    assert official_path.num_steps == 2
    assert official_path.depth == 2
    assert [node.smiles for node in official_path.get_all_steps()] == ["CCO", "CC=O"]
    assert {node.smiles for node in official_path.get_starting_materials()} == {
        "CC",
        "O",
    }
    assert PathConsistencyChecker().check(official_path) is True
    assert [
        variant.num_steps for variant in PathCollapser().collapse(official_path)
    ] == [
        2,
        1,
    ]


@pytest.mark.parametrize(
    ("pathway", "message"),
    [
        ([], "at least one step"),
        ([{"product": "not-smiles", "reactants": ["CC"]}], "Invalid product"),
        ([{"product": "CCO", "reactants": ["not-smiles"]}], "Invalid reactant"),
        (
            [
                {"product": "CCO", "reactants": ["CC"]},
                {"product": "C(C)O", "reactants": ["C"]},
            ],
            "Duplicate producing steps",
        ),
        (
            [
                {"product": "CCO", "reactants": ["CC"]},
                {"product": "CCN", "reactants": ["CN"]},
            ],
            "exactly one root",
        ),
        (
            [
                {"product": "CCO", "reactants": ["CCN"]},
                {"product": "CCN", "reactants": ["CCO"]},
            ],
            "cycle",
        ),
        (
            [
                {"product": "CCO", "reactants": ["CC"]},
                {"product": "CCN", "reactants": ["CNC"]},
                {"product": "CNC", "reactants": ["CCN"]},
            ],
            "disconnected or cyclic",
        ),
        (
            [
                {"product": "CCOC", "reactants": ["CCO", "CO"]},
                {"product": "CCO", "reactants": ["CC=O"]},
                {"product": "CO", "reactants": ["CCO", "O"]},
            ],
            "consumed more than once",
        ),
        (
            [
                {
                    "products": [{"smiles": "CCO"}, {"smiles": "CCN"}],
                    "reactants": [{"smiles": "CC"}],
                }
            ],
            "exactly one product",
        ),
        (
            [
                {
                    "product": "CCO",
                    "product_smiles": "CCN",
                    "reactants": ["CC"],
                }
            ],
            "ambiguous product fields",
        ),
        (
            [
                {
                    "product": "CCO",
                    "reactants": ["CC"],
                    "reactants_smiles": ["CN"],
                }
            ],
            "ambiguous reactant fields",
        ),
        ([{"product": " ", "reactants": ["CC"]}], "Invalid product"),
        ([{"product": "CCO", "reactants": [" "]}], "missing reactant"),
    ],
)
def test_invalid_route_shapes_are_rejected(
    pathway: list[dict[str, Any]],
    message: str,
) -> None:
    with pytest.raises(UrsaRouteError, match=message):
        _to_retrocast_route(pathway)


def test_invalid_sqlite_is_reported_by_the_real_backend(tmp_path: Path) -> None:
    config = existing_asset_config(tmp_path)
    config.chemcensor_db_path.write_bytes(b"not a sqlite database")

    with pytest.raises(UrsaConfigurationError, match="initialize URSA") as exc_info:
        UrsaAdapter(config).annotate(LINEAR_PATHWAY)

    assert exc_info.value.__cause__ is not None


@pytest.mark.slow
def test_score_pathway_with_real_ursa_when_official_assets_are_configured() -> None:
    result = score.score_pathway(
        RELEASED_EXAMPLE_PATHWAY,
        ursa_config=configured_real_assets(),
    )
    annotation = result["ursa"]

    assert result["route_summary"]["n_steps"] == 7
    assert annotation["path_id"] == "deepretro-route"
    assert annotation["is_consistent"] is True
    assert annotation["all_bb_found"] is True
    assert annotation["passes_solv_0"] is True
    assert annotation["passes_solv_1"] is True
    assert annotation["passes_solv_2"] is True

    for key in ("best_variant_solv_1", "best_variant_solv_2"):
        steps = annotation[key]["steps"]
        assert steps
        assert RELEASED_EXAMPLE_ROOT_REACTION in {
            step["reaction_smiles"] for step in steps
        }
        assert all(isinstance(step["score_without_fg"], (int, float)) for step in steps)
        assert all(isinstance(step["score_with_fg"], (int, float)) for step in steps)
