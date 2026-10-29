"""Shared pytest fixtures for DeepRetro unit tests."""

from __future__ import annotations

import importlib.util
import os
import socket
import sys
from collections.abc import Callable
from pathlib import Path
from typing import NoReturn, Optional
from types import ModuleType

import pytest

import deepretro

from deepretro.metadata_types import (
    ConditionsPayload,
    ConditionsRecommender,
    ConditionsStatusPayload,
    LiteratureRecommender,
    LiteratureStatusPayload,
    MoleculeRecord,
    ReagentRecommender,
    ReagentStatusPayload,
)
from deepretro.utils.parse import RetrosynthesisRouteParser
from deepretro.utils.utils_molecule import calc_chemical_formula, calc_mol_wt

AZ_MODULE_PATH = Path(deepretro.__file__).resolve().parent / "utils" / "az.py"
AZ_MODULE_NAME = "_deepretro_utils_az_under_test"


# ---------------------------------------------------------------------------
# AiZynthFinder fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def az_models_dir() -> Path:
    """Use explicitly provisioned AiZynthFinder models for integration tests.

    Set ``AZ_TEST_MODELS_DIR`` to a directory containing ``config.yml``.
    Tests never download model files themselves.
    """
    env_dir = os.environ.get("AZ_TEST_MODELS_DIR")
    if not env_dir:
        pytest.skip("Set AZ_TEST_MODELS_DIR to run AiZynthFinder integration tests")
    models_dir = Path(env_dir).expanduser().resolve()
    if not (models_dir / "config.yml").is_file():
        pytest.fail(f"AZ_TEST_MODELS_DIR has no config.yml: {models_dir}")
    pytest.importorskip("aizynthfinder")
    return models_dir


@pytest.fixture
def az_module() -> ModuleType:
    """Load the az module from source."""
    sys.modules.pop(AZ_MODULE_NAME, None)
    spec = importlib.util.spec_from_file_location(AZ_MODULE_NAME, AZ_MODULE_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[AZ_MODULE_NAME] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def az_module_with_models(az_module: ModuleType, az_models_dir: Path) -> ModuleType:
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


# ---------------------------------------------------------------------------
# Metadata recommender test doubles
# ---------------------------------------------------------------------------


@pytest.fixture
def recommender_calls() -> list[str]:
    """Shared call-tracking list for recommender spy fixtures."""
    return []


@pytest.fixture
def spy_reagent_recommender(
    recommender_calls: list[str],
) -> ReagentRecommender:
    """Reagent recommender that logs calls and returns a water reagent."""

    def _recommender(
        reactants: list[MoleculeRecord],
        product: list[MoleculeRecord],
        model: str,
        temperature: float,
    ) -> ReagentStatusPayload:
        recommender_calls.append(f"reagent:{model}:{temperature}")
        return 200, [{"smiles": "O", "reagent_metadata": {"name": ""}}]

    return _recommender


@pytest.fixture
def spy_conditions_recommender(
    recommender_calls: list[str],
) -> ConditionsRecommender:
    """Conditions recommender that logs calls and returns standard conditions."""

    def _recommender(
        reactants: list[MoleculeRecord],
        product: list[MoleculeRecord],
        reagents: list[MoleculeRecord],
        model: str,
        temperature: float,
    ) -> ConditionsStatusPayload:
        recommender_calls.append(f"conditions:{model}:{temperature}")
        return 200, {
            "temperature": "25 C",
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 h",
        }

    return _recommender


@pytest.fixture
def spy_literature_recommender(
    recommender_calls: list[str],
) -> LiteratureRecommender:
    """Literature recommender that logs calls and returns a stub DOI."""

    def _recommender(
        reactants: list[MoleculeRecord],
        product: list[MoleculeRecord],
        reagents: list[MoleculeRecord],
        conditions: ConditionsPayload,
        model: str,
        temperature: float,
    ) -> LiteratureStatusPayload:
        recommender_calls.append(f"literature:{model}:{temperature}")
        return 200, {"doi": "10.1000/example"}

    return _recommender


@pytest.fixture
def failing_reagent_recommender() -> ReagentRecommender:
    """Reagent recommender that always returns a 404 failure."""

    def _recommender(
        reactants: list[MoleculeRecord],
        product: list[MoleculeRecord],
        model: str,
        temperature: float,
    ) -> ReagentStatusPayload:
        return 404, ""

    return _recommender


@pytest.fixture(autouse=True)
def _offline_langfuse_events(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep the suite offline: never build a real Langfuse client in tests.

    Tests that want to observe events monkeypatch ``_langfuse_client`` again
    with a recorder.
    """
    from deepretro.utils import llm_trace

    monkeypatch.setattr(llm_trace, "_langfuse_client", lambda: None)


def pytest_addoption(parser: pytest.Parser) -> None:
    """Register the explicit opt-in for network and model integration tests."""
    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="Run integration tests that may call paid providers or use local models",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Skip integration tests unless the caller explicitly enables them."""
    if config.getoption("--run-integration"):
        return
    skip = pytest.mark.skip(reason="Integration test: pass --run-integration to enable")
    for item in items:
        if item.get_closest_marker("integration") is not None:
            item.add_marker(skip)


@pytest.fixture(autouse=True)
def _block_unit_test_network(
    request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Fail unexpected socket connections in unit tests, even if keys exist.

    This guard covers Python sockets in the test process. External programs
    and subprocesses must still use explicit test doubles.
    """
    if request.node.get_closest_marker("integration") is not None:
        return

    def blocked(*args: object, **kwargs: object) -> NoReturn:
        raise RuntimeError(
            "Network access is disabled in unit tests. Use a test double or mark "
            "the test integration and pass --run-integration."
        )

    monkeypatch.setattr(socket, "create_connection", blocked)
    monkeypatch.setattr(socket.socket, "connect", blocked)
    monkeypatch.setattr(socket.socket, "connect_ex", blocked)
