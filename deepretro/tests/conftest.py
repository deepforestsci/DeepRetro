"""Shared pytest fixtures for DeepRetro unit tests."""

from __future__ import annotations

from collections.abc import Callable
from typing import NoReturn, Optional

import pytest

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
