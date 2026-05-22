"""Shared pytest fixtures for DeepRetro unit tests."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import NoReturn, Optional

import pytest

from deepretro.utils.parse import RetrosynthesisRouteParser
from deepretro.utils.utils_molecule import calc_chemical_formula, calc_mol_wt

PROJECT_ROOT = Path(__file__).resolve().parents[2]
AZ_MODULE_PATH = PROJECT_ROOT / "deepretro" / "utils" / "az.py"
AZ_MODULE_NAME = "_deepretro_utils_az_under_test"


# ---------------------------------------------------------------------------
# AiZynthFinder fixtures
# ---------------------------------------------------------------------------


def _download_aizynth_models(models_dir: Path) -> None:
    """Download AiZynthFinder public models to *models_dir* (skips if config.yml exists)."""
    if (models_dir / "config.yml").exists():
        return
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
        timeout=600,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to download AiZynthFinder models: {result.stderr or result.stdout}"
        )


@pytest.fixture(scope="session")
def az_models_dir(tmp_path_factory):
    """Session-scoped fixture providing a directory with AiZynthFinder models.

    Honours ``AZ_TEST_MODELS_DIR`` so CI can point to a cached directory
    and skip the download entirely.
    """
    env_dir = os.environ.get("AZ_TEST_MODELS_DIR")
    if env_dir:
        models_dir = Path(env_dir)
    else:
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


# ---------------------------------------------------------------------------
# Reaction classifier / parse fixtures
# ---------------------------------------------------------------------------


class DummyReactionClassifier:
    """Minimal classifier stub returning one reaction encoding."""

    def predict(self, fingerprints: list[list[int]]) -> list[int]:
        """Return a single deterministic reaction encoding."""
        assert fingerprints == [[1, 0, 1, 0]]
        return [0]


class UnknownReactionClassifier:
    """Classifier stub returning an unmapped reaction encoding."""

    def predict(self, fingerprints: list[list[int]]) -> list[int]:
        """Return an encoding that is intentionally absent from the mapping."""
        return [999]


@pytest.fixture
def dummy_reaction_classifier() -> DummyReactionClassifier:
    """Return a deterministic classifier test double."""
    return DummyReactionClassifier()


@pytest.fixture
def dummy_model_loader(
    dummy_reaction_classifier: DummyReactionClassifier,
) -> Callable[[str], DummyReactionClassifier]:
    """Return a model loader that reuses the dummy classifier."""

    def _loader(path: str) -> DummyReactionClassifier:
        return dummy_reaction_classifier

    return _loader


@pytest.fixture
def unknown_model_loader() -> Callable[[str], UnknownReactionClassifier]:
    """Return a model loader for an unmapped classifier output."""

    def _loader(path: str) -> UnknownReactionClassifier:
        return UnknownReactionClassifier()

    return _loader


@pytest.fixture
def missing_model_loader() -> Callable[[str], NoReturn]:
    """Return a model loader that simulates an absent model file."""

    def _loader(path: str) -> NoReturn:
        raise FileNotFoundError

    return _loader


@pytest.fixture
def missing_model_loader_with_message() -> Callable[[str], NoReturn]:
    """Return a missing-model loader with a stable error message."""

    def _loader(path: str) -> NoReturn:
        raise FileNotFoundError("model.joblib not found")

    return _loader


@pytest.fixture
def fingerprint_calculator() -> Callable[[str], list[int]]:
    """Return a deterministic fingerprint calculator."""

    def _calculator(smiles: str) -> list[int]:
        return [1, 0]

    return _calculator


@pytest.fixture
def missing_fingerprint_calculator() -> Callable[[str], Optional[list[int]]]:
    """Return a fingerprint calculator that simulates RDKit failures."""

    def _calculator(smiles: str) -> Optional[list[int]]:
        return None

    return _calculator


@pytest.fixture
def route_parser_factory() -> Callable[..., RetrosynthesisRouteParser]:
    """Return a factory for parsers with deterministic defaults."""

    def _factory(
        scalability_calculator: Optional[Callable[[str, str], int]] = None,
        basic_molecules: Optional[set[str]] = None,
    ) -> RetrosynthesisRouteParser:
        return RetrosynthesisRouteParser(
            basic_molecules=set() if basic_molecules is None else basic_molecules,
            chemical_formula_calculator=calc_chemical_formula,
            mass_calculator=calc_mol_wt,
            scalability_calculator=(
                (lambda reactant, product: 0)
                if scalability_calculator is None
                else scalability_calculator
            ),
        )

    return _factory
