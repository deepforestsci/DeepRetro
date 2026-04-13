"""Unit tests for deepretro/utils/az.py.

These tests run only when aizynthfinder is installed. They use real AiZynthFinder
with models downloaded to a temporary directory at session start. No mocking.

Speed optimizations:
- Tests that run full inference are marked @pytest.mark.slow; skip with -m "not slow"
"""

from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest
from deepretro.utils.cache import CacheManager

PROJECT_ROOT = Path(__file__).resolve().parents[2]
AZ_MODULE_PATH = PROJECT_ROOT / "deepretro" / "utils" / "az.py"
AZ_MODULE_NAME = "_deepretro_utils_az_under_test"
_BASIC_ROUTE = [
    {
        "type": "mol",
        "hide": False,
        "smiles": "CCO",
        "is_chemical": True,
        "in_stock": True,
    }
]


def _download_aizynth_models(models_dir: Path) -> None:
    """Download AiZynthFinder public models to the given directory."""
    models_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "aizynthfinder.tools.download_public_data",
            str(models_dir),
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to download AiZynthFinder models: {result.stderr or result.stdout}"
        )


@pytest.fixture(scope="session")
def az_models_dir(tmp_path_factory):
    """Session-scoped fixture: download AiZynthFinder models to a temp directory."""
    pytest.importorskip("aizynthfinder")
    models_dir = tmp_path_factory.mktemp("aizynth_models")
    _download_aizynth_models(models_dir)
    return models_dir


@pytest.fixture
def az_module():
    """Load the az module from source."""
    sys.modules.pop(AZ_MODULE_NAME, None)
    spec = importlib.util.spec_from_file_location(AZ_MODULE_NAME, AZ_MODULE_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[AZ_MODULE_NAME] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def az_module_with_models(az_module, az_models_dir):
    """Az module with paths configured to use downloaded models."""
    config_path = az_models_dir / "config.yml"
    az_module.AZ_MODEL_CONFIG_PATH = str(config_path)
    az_module.AZ_MODELS_PATH = str(az_models_dir)
    return az_module


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
def test_is_basic_molecule_chemical_thresholds(az_module, smiles, expected):
    """Test is_basic_molecule chemical thresholds."""
    assert az_module.is_basic_molecule(smiles) is expected


@pytest.mark.slow
def test_run_az_returns_valid_result(az_module_with_models):
    """Test run_az with real AiZynthFinder returns (bool, list) with valid structure."""
    status, routes = az_module_with_models.run_az("C1CCCCC1", az_model="USPTO")
    assert isinstance(status, bool)
    assert isinstance(routes, list)
    # Routes may be empty or contain dicts with expected keys
    for route in routes:
        assert isinstance(route, dict)
        assert "type" in route or "smiles" in route


@pytest.mark.slow
def test_run_az_uses_fallback_config_when_model_missing(az_module_with_models):
    """Test run_az falls back to AZ_MODEL_CONFIG_PATH when model-specific config missing."""
    # MISSING_MODEL has no config in AZ_MODELS_PATH, so uses fallback
    status, routes = az_module_with_models.run_az("C1CCCCC1", az_model="MISSING_MODEL")
    assert isinstance(status, bool)
    assert isinstance(routes, list)


def test_run_az_raises_if_no_config_available(tmp_path, monkeypatch, az_module):
    """Test run_az raises FileNotFoundError when no config exists."""
    missing_models_root = tmp_path / "models"
    missing_fallback = tmp_path / "missing.yml"
    monkeypatch.setattr(az_module, "AZ_MODELS_PATH", str(missing_models_root))
    monkeypatch.setattr(az_module, "AZ_MODEL_CONFIG_PATH", str(missing_fallback))
    monkeypatch.setattr(az_module, "BASIC_MOLECULES", [])

    with pytest.raises(FileNotFoundError, match=re.escape(str(missing_fallback))):
        az_module.run_az("CCCCC", az_model="MISSING_MODEL")


def test_run_az_short_circuits_for_feedstock_smiles(az_module):
    """Test run_az returns early for basic molecules without calling AiZynthFinder."""
    status, routes = az_module.run_az("CCO", az_model="USPTO")
    assert status is True
    assert routes == _BASIC_ROUTE


def test_run_az_has_no_cache_side_effect_without_cache(az_module):
    """run_az should not create cache state unless the caller passes one in."""
    cache = CacheManager()

    status, routes = az_module.run_az("CCO", az_model="USPTO")

    assert status is True
    assert routes == _BASIC_ROUTE
    assert cache.stats().num_entries == 0


def test_run_az_uses_explicit_cache(az_module):
    """run_az should reuse results only when an explicit cache is provided."""
    cache = CacheManager()

    first = az_module.run_az("CCO", az_model="USPTO", cache=cache)
    second = az_module.run_az("CCO", az_model="USPTO", cache=cache)

    assert first == second
    assert cache.stats().hits == 1
    assert cache.stats().misses == 1
    assert cache.stats().num_entries == 1
    cache.close()


@pytest.mark.slow
def test_run_az_with_img_returns_valid_result(az_module_with_models):
    """Test run_az_with_img with real AiZynthFinder returns (bool, list, images)."""
    status, routes, images = az_module_with_models.run_az_with_img("C1CCCCC1")
    assert isinstance(status, bool)
    assert isinstance(routes, list)
    # images can be list of PIL Images or None
    if images is not None:
        assert isinstance(images, (list, tuple))


def test_run_az_with_img_short_circuits_for_basic_molecules(az_module):
    """Test run_az_with_img returns early for basic molecules."""
    status, routes, images = az_module.run_az_with_img("CCO")
    assert status is True
    assert routes == _BASIC_ROUTE
    assert images is None


def test_run_az_with_img_uses_explicit_cache(az_module):
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
