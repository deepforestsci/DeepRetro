"""Unit tests for deepretro/utils/az.py.

Most tests use feedstock shortcuts and need no optional backend. Tests marked
``integration`` and ``slow`` use real models configured by ``AZ_TEST_MODELS_DIR``.
They require the explicit ``--run-integration`` pytest option.
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType

import re

import pytest
from deepretro.utils.cache import CacheManager

_BASIC_ROUTE = [
    {
        "type": "mol",
        "hide": False,
        "smiles": "CCO",
        "is_chemical": True,
        "in_stock": True,
    }
]


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("C", True),
        ("O=S(=O)(O)O", True),
        ("CCCCC", False),
        ("c1ccccc1", False),
        ("not_a_smiles!!!", False),
    ],
)
def test_is_basic_molecule_chemical_thresholds(
    az_module: ModuleType, smiles: str, expected: bool
) -> None:
    """Test is_basic_molecule chemical thresholds."""
    assert az_module.is_basic_molecule(smiles) is expected


# ---------------------------------------------------------------------------
# _run_az_core tests
# ---------------------------------------------------------------------------


def test_run_az_core_short_circuits_basic_molecule(az_module: ModuleType) -> None:
    """_run_az_core returns (True, route, None) for basic molecules."""
    status, routes, finder = az_module._run_az_core("CCO")
    assert status is True
    assert routes == _BASIC_ROUTE
    assert finder is None


def test_run_az_core_short_circuits_with_az_model(az_module: ModuleType) -> None:
    """az_model param is irrelevant for basic molecules."""
    status, routes, finder = az_module._run_az_core("CCO", az_model="USPTO")
    assert status is True
    assert routes == _BASIC_ROUTE
    assert finder is None


def test_run_az_core_raises_file_not_found_before_import_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, az_module: ModuleType
) -> None:
    """Config resolution should fail before the AiZynthFinder import check."""
    missing_models_root = tmp_path / "models"
    missing_fallback = tmp_path / "missing.yml"
    monkeypatch.setattr(az_module, "AZ_MODELS_PATH", str(missing_models_root))
    monkeypatch.setattr(az_module, "AZ_MODEL_CONFIG_PATH", str(missing_fallback))
    monkeypatch.setattr(az_module, "BASIC_MOLECULES", [])

    with pytest.raises(FileNotFoundError, match=re.escape(str(missing_fallback))):
        az_module._run_az_core("CCCCC", az_model="MISSING_MODEL")


@pytest.mark.integration
@pytest.mark.slow
def test_run_az_core_returns_finder_for_real_molecule(
    az_module_with_models: ModuleType,
) -> None:
    """_run_az_core returns a non-None finder for a real molecule."""
    status, routes, finder = az_module_with_models._run_az_core(
        "CC(=O)Oc1ccccc1C(=O)O", az_model="USPTO"
    )
    assert isinstance(status, bool)
    assert isinstance(routes, list)
    assert finder is not None


@pytest.mark.integration
@pytest.mark.slow
def test_run_az_returns_valid_result(az_module_with_models: ModuleType) -> None:
    """Test run_az with real AiZynthFinder returns (bool, list) with valid structure."""
    status, routes = az_module_with_models.run_az("C1CCCCC1", az_model="USPTO")
    assert isinstance(status, bool)
    assert isinstance(routes, list)
    # Routes may be empty or contain dicts with expected keys
    for route in routes:
        assert isinstance(route, dict)
        assert "type" in route or "smiles" in route


@pytest.mark.integration
@pytest.mark.slow
def test_run_az_uses_fallback_config_when_model_missing(
    az_module_with_models: ModuleType,
) -> None:
    """Test run_az falls back to AZ_MODEL_CONFIG_PATH when model-specific config missing."""
    # MISSING_MODEL has no config in AZ_MODELS_PATH, so uses fallback
    status, routes = az_module_with_models.run_az("C1CCCCC1", az_model="MISSING_MODEL")
    assert isinstance(status, bool)
    assert isinstance(routes, list)


def test_run_az_raises_if_no_config_available(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, az_module: ModuleType
) -> None:
    """Test run_az raises FileNotFoundError when no config exists."""
    missing_models_root = tmp_path / "models"
    missing_fallback = tmp_path / "missing.yml"
    monkeypatch.setattr(az_module, "AZ_MODELS_PATH", str(missing_models_root))
    monkeypatch.setattr(az_module, "AZ_MODEL_CONFIG_PATH", str(missing_fallback))
    monkeypatch.setattr(az_module, "BASIC_MOLECULES", [])

    with pytest.raises(FileNotFoundError, match=re.escape(str(missing_fallback))):
        az_module.run_az("CCCCC", az_model="MISSING_MODEL")


def test_run_az_short_circuits_for_feedstock_smiles(az_module: ModuleType) -> None:
    """Test run_az returns early for basic molecules without calling AiZynthFinder."""
    status, routes = az_module.run_az("CCO", az_model="USPTO")
    assert status is True
    assert routes == _BASIC_ROUTE


def test_run_az_has_no_cache_side_effect_without_cache(az_module: ModuleType) -> None:
    """run_az should not create cache state unless the caller passes one in."""
    cache = CacheManager()

    status, routes = az_module.run_az("CCO", az_model="USPTO")

    assert status is True
    assert routes == _BASIC_ROUTE
    assert cache.stats().num_entries == 0


def test_run_az_uses_explicit_cache(az_module: ModuleType) -> None:
    """run_az should reuse results only when an explicit cache is provided."""
    cache = CacheManager()

    first = az_module.run_az("CCO", az_model="USPTO", cache=cache)
    second = az_module.run_az("CCO", az_model="USPTO", cache=cache)

    assert first == second
    assert cache.stats().hits == 1
    assert cache.stats().misses == 1
    assert cache.stats().num_entries == 1
    cache.close()


@pytest.mark.integration
@pytest.mark.slow
def test_run_az_with_img_returns_valid_result(
    az_module_with_models: ModuleType,
) -> None:
    """Test run_az_with_img with real AiZynthFinder returns (bool, list, images)."""
    status, routes, images = az_module_with_models.run_az_with_img("C1CCCCC1")
    assert isinstance(status, bool)
    assert isinstance(routes, list)
    # images can be list of PIL Images or None
    if images is not None:
        assert isinstance(images, (list, tuple))


def test_run_az_with_img_short_circuits_for_basic_molecules(
    az_module: ModuleType,
) -> None:
    """Test run_az_with_img returns early for basic molecules."""
    status, routes, images = az_module.run_az_with_img("CCO")
    assert status is True
    assert routes == _BASIC_ROUTE
    assert images is None


def test_run_az_with_img_uses_explicit_cache(az_module: ModuleType) -> None:
    """run_az_with_img should also cache only through the provided instance."""
    cache = CacheManager()

    first = az_module.run_az_with_img("CCO", cache=cache)
    second = az_module.run_az_with_img("CCO", cache=cache)

    assert first == second
    assert first == (True, _BASIC_ROUTE, None)
    assert cache.stats().hits == 1
    assert cache.stats().misses == 1
    assert cache.stats().num_entries == 1
    cache.close()


@pytest.mark.parametrize("variable", ["AZ_MODEL_CONFIG_PATH", "AZ_MODELS_PATH"])
def test_absolute_environment_path_is_preserved(
    variable: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    az_module: ModuleType,
) -> None:
    """Absolute model paths must not be prefixed by the checkout directory."""
    import runpy

    path = tmp_path / "models"
    config = path / "USPTO" / "config.yml" if variable == "AZ_MODELS_PATH" else path
    config.parent.mkdir(parents=True, exist_ok=True)
    config.write_text("{}", encoding="utf-8")
    monkeypatch.setenv(variable, str(path))
    namespace = runpy.run_path(az_module.__file__)
    assert namespace[variable] == str(path)
    model = "USPTO" if variable == "AZ_MODELS_PATH" else None
    assert namespace["_resolve_config"](model) == str(config)


def test_relative_environment_path_is_checkout_relative(
    monkeypatch: pytest.MonkeyPatch, az_module: ModuleType
) -> None:
    """Existing relative-path configurations keep their source-checkout meaning."""
    monkeypatch.setenv("AZ_MODELS_PATH", "models")
    assert az_module._configured_path("AZ_MODELS_PATH") == str(
        az_module.PROJECT_ROOT / "models"
    )


def test_unset_environment_path_is_empty(
    monkeypatch: pytest.MonkeyPatch, az_module: ModuleType
) -> None:
    """Missing configuration must not become a literal path ending in None."""
    monkeypatch.delenv("AZ_MODELS_PATH", raising=False)
    assert az_module._configured_path("AZ_MODELS_PATH") == ""


def test_single_step_ranks_templates_and_skips_failed_applications(
    monkeypatch: pytest.MonkeyPatch, az_module: ModuleType
) -> None:
    """A failed high-prior template must not hide later usable suggestions."""
    import sys
    from types import SimpleNamespace
    from typing import Any

    class FailedAction:
        @property
        def reactants(self) -> list[Any]:
            raise ValueError("template does not apply")

    chem = ModuleType("aizynthfinder.chem")
    chem.TreeMolecule = lambda **kwargs: SimpleNamespace(**kwargs)
    monkeypatch.setitem(sys.modules, "aizynthfinder.chem", chem)
    actions = [
        SimpleNamespace(
            reactants=[[SimpleNamespace(smiles="CC"), SimpleNamespace(smiles="O")]]
        ),
        FailedAction(),
        SimpleNamespace(reactants=[]),
    ]
    finder = SimpleNamespace(
        expansion_policy=SimpleNamespace(
            get_actions=lambda molecules: (actions, [0.5, 0.9, 0.2])
        )
    )
    monkeypatch.setattr(az_module, "_resolve_config", lambda model: "example.yml")
    monkeypatch.setattr(az_module, "_get_finder", lambda path: finder)

    assert az_module.az_single_step("CCO", top_k=3) == [
        {"precursors": ["CC", "O"], "probability": 0.5}
    ]
    assert az_module.az_single_step("CCO", top_k=1) == []
