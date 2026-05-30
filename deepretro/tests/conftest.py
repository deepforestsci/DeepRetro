"""Shared pytest fixtures for DeepRetro unit tests."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, NoReturn, Optional

import pytest

from deepretro.algorithms.autosolve import reaction_tree, unsolved_leaf
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


# ---------------------------------------------------------------------------
# AutoSolver test doubles
# ---------------------------------------------------------------------------

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SALICYLIC_ACID = "OC(=O)c1ccccc1O"
ACETIC_ANHYDRIDE = "CC(=O)OC(C)=O"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
PARA_AMINOPHENOL = "Nc1ccc(O)cc1"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
ACETIC_ACID = "CC(=O)O"


def solved_route(smiles: str) -> dict[str, Any]:
    """Return a minimal solved route node."""
    return {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": True,
    }


def az_always_fails(smiles: str, az_model: str) -> tuple[bool, list[Any]]:
    """AZ runner that never solves any molecule."""
    return False, []


def llm_returns_nothing(molecule: str, **kwargs: Any) -> tuple[list, list, list]:
    """LLM runner that returns no pathways."""
    return [], [], []


def aspirin_hydrolysis_tree() -> dict[str, Any]:
    """Aspirin -> salicylic acid + acetic anhydride (one reaction step)."""
    return reaction_tree(
        ASPIRIN,
        [unsolved_leaf(SALICYLIC_ACID), unsolved_leaf(ACETIC_ANHYDRIDE)],
        [0.9],
    )
