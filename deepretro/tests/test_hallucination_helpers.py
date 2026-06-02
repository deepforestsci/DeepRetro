"""Tests for hallucination checker adapter helpers."""

from __future__ import annotations


import pytest

from deepretro.models.hallucination_helpers import (
    build_ml_checker,
    filter_with_checker,
    resolve_hallucination,
)

# Real molecule SMILES for meaningful tests
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SALICYLIC_ACID = "OC(=O)c1ccccc1O"
ACETIC_ACID = "CC(=O)O"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
PARA_AMINOPHENOL = "Nc1ccc(O)cc1"


# ---------------------------------------------------------------------------
# Classifier test doubles — injected, not mocked
# ---------------------------------------------------------------------------


class _AlwaysCleanClassifier:
    """Classifier that marks every reaction as non-hallucinated."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []

    def predict_single(self, product: str, reactants: str) -> dict[str, bool]:
        self.calls.append((product, reactants))
        return {"is_hallucination": False}


class _SelectiveClassifier:
    """Classifier that marks specific reactants as hallucinated."""

    def __init__(self, hallucinated_reactants: set[str]) -> None:
        self.hallucinated = hallucinated_reactants
        self.calls: list[tuple[str, str]] = []

    def predict_single(self, product: str, reactants: str) -> dict[str, bool]:
        self.calls.append((product, reactants))
        return {"is_hallucination": reactants in self.hallucinated}


# ---------------------------------------------------------------------------
# build_ml_checker tests
# ---------------------------------------------------------------------------


class TestBuildMlChecker:
    def test_keeps_non_hallucinated_aspirin_precursors(self) -> None:
        """Both real precursors pass the classifier — both retained."""
        classifier = _AlwaysCleanClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID], [ACETIC_ACID]])

        assert status == 200
        assert retained == [[SALICYLIC_ACID], [ACETIC_ACID]]
        assert len(classifier.calls) == 2
        assert classifier.calls[0] == (ASPIRIN, SALICYLIC_ACID)
        assert classifier.calls[1] == (ASPIRIN, ACETIC_ACID)

    def test_filters_hallucinated_pathway(self) -> None:
        """Classifier flags one precursor — only the clean one survives."""
        classifier = _SelectiveClassifier(hallucinated_reactants={SALICYLIC_ACID})
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID], [ACETIC_ACID]])

        assert status == 200
        assert retained == [[ACETIC_ACID]]

    def test_drops_invalid_smiles_before_classifier(self) -> None:
        """Invalid SMILES are filtered out before reaching the classifier."""
        classifier = _AlwaysCleanClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [["not>>>valid"]])

        assert status == 200
        assert retained == []
        assert classifier.calls == []

    def test_missing_hallucination_key_defaults_to_filtered(self) -> None:
        """Classifier returning no is_hallucination key treats pathway as hallucinated."""

        class _NoKeyClassifier:
            def predict_single(self, product: str, reactants: str) -> dict:
                return {"confidence": 0.5}

        checker = build_ml_checker(_NoKeyClassifier())
        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID]])

        assert status == 200
        assert retained == []

    def test_multi_reactant_pathway_joins_with_dot(self) -> None:
        """Multiple reactants in one pathway are joined with '.' for the classifier."""
        classifier = _AlwaysCleanClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(PARACETAMOL, [[PARA_AMINOPHENOL, ACETIC_ACID]])

        assert status == 200
        assert retained == [[PARA_AMINOPHENOL, ACETIC_ACID]]
        expected_joined = f"{PARA_AMINOPHENOL}.{ACETIC_ACID}"
        assert classifier.calls == [(PARACETAMOL, expected_joined)]


# ---------------------------------------------------------------------------
# filter_with_checker tests
# ---------------------------------------------------------------------------


class TestFilterWithChecker:
    def test_preserves_alignment_after_filtering(self) -> None:
        """Retained pathway keeps its original explanation and confidence."""

        def keeps_only_second(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 200, [pathways[1]]

        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN,
            [[SALICYLIC_ACID], [ACETIC_ACID]],
            ["hydrolysis", "deacetylation"],
            [0.3, 0.9],
            keeps_only_second,
        )

        assert pathways == [[ACETIC_ACID]]
        assert explanations == ["deacetylation"]
        assert confidence == [0.9]

    def test_no_checker_passes_through(self) -> None:
        """None checker returns everything unchanged."""
        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN,
            [[SALICYLIC_ACID], [ACETIC_ACID]],
            ["first", "second"],
            [0.1, 0.9],
            None,
        )

        assert pathways == [[SALICYLIC_ACID], [ACETIC_ACID]]
        assert explanations == ["first", "second"]
        assert confidence == [0.1, 0.9]

    def test_non_200_status_drops_everything(self) -> None:
        """Non-200 checker status discards all pathways and metadata."""

        def fails(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 500, pathways

        result = filter_with_checker(
            ASPIRIN, [[SALICYLIC_ACID]], ["only"], [0.8], fails
        )

        assert result == ([], [], [])

    def test_duplicate_pathways_consume_metadata_in_order(self) -> None:
        """When checker returns duplicates, metadata is consumed sequentially."""

        def returns_duplicates(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 200, [[SALICYLIC_ACID], [SALICYLIC_ACID]]

        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN,
            [[SALICYLIC_ACID], [SALICYLIC_ACID], [ACETIC_ACID]],
            ["first", "second", "third"],
            [0.1, 0.2, 0.3],
            returns_duplicates,
        )

        assert pathways == [[SALICYLIC_ACID], [SALICYLIC_ACID]]
        assert explanations == ["first", "second"]
        assert confidence == [0.1, 0.2]

    def test_empty_pathways_short_circuits(self) -> None:
        """Empty pathway list skips checker entirely."""
        call_count = 0

        def should_not_be_called(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            nonlocal call_count
            call_count += 1
            return 200, pathways

        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN, [], [], [], should_not_be_called
        )

        assert pathways == []
        assert call_count == 0

    def test_unknown_pathway_from_checker_is_dropped(self) -> None:
        """Checker returning a pathway not in the original input is silently dropped."""

        def injects_unknown(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 200, [[PARACETAMOL]]  # not in original input

        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN,
            [[SALICYLIC_ACID]],
            ["hydrolysis"],
            [0.9],
            injects_unknown,
        )

        assert pathways == []
        assert explanations == []
        assert confidence == []


# ---------------------------------------------------------------------------
# resolve_hallucination tests
# ---------------------------------------------------------------------------


class TestResolveHallucination:
    def test_none_mode_returns_none(self) -> None:
        assert resolve_hallucination("none", None) is None

    def test_heuristic_mode_returns_none(self) -> None:
        assert resolve_hallucination("heuristic", None) is None

    def test_ml_mode_with_classifier_returns_callable(self) -> None:
        classifier = _AlwaysCleanClassifier()
        checker = resolve_hallucination("ml", classifier)
        assert callable(checker)

    def test_ml_mode_without_classifier_raises(self) -> None:
        with pytest.raises(ValueError, match="requires"):
            resolve_hallucination("ml", None)

    def test_unknown_mode_raises(self) -> None:
        with pytest.raises(ValueError, match="hallucination_mode"):
            resolve_hallucination("unknown", None)

    def test_case_insensitive(self) -> None:
        assert resolve_hallucination("NONE", None) is None
        assert resolve_hallucination("Heuristic", None) is None

    def test_ml_checker_actually_filters(self) -> None:
        """End-to-end: resolve → build → call with real molecules."""
        classifier = _SelectiveClassifier(hallucinated_reactants={SALICYLIC_ACID})
        checker = resolve_hallucination("ml", classifier)
        assert checker is not None

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID], [ACETIC_ACID]])

        assert status == 200
        assert retained == [[ACETIC_ACID]]
