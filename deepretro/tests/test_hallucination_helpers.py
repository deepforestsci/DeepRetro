"""Tests for hallucination checker adapter helpers."""

from __future__ import annotations

from typing import Any

import pytest

from deepretro.models.hallucination_helpers import (
    build_ml_checker,
    filter_with_checker,
    resolve_hallucination,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SALICYLIC_ACID = "OC(=O)c1ccccc1O"
ACETIC_ACID = "CC(=O)O"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
PARA_AMINOPHENOL = "Nc1ccc(O)cc1"


class StubClassifier:
    """Configurable classifier double that records calls."""

    def __init__(self, responses: list[dict[str, Any]] | None = None) -> None:
        self.responses = list(responses) if responses else []
        self.calls: list[tuple[str, str]] = []

    def predict_single(
        self, product_smiles: str, reactants_smiles: str
    ) -> dict[str, Any]:
        self.calls.append((product_smiles, reactants_smiles))
        if self.responses:
            return self.responses.pop(0)
        return {"is_hallucination": False}


class SelectiveClassifier:
    """Classifier that flags specific reactant SMILES as hallucinated."""

    def __init__(self, hallucinated: set[str]) -> None:
        self.hallucinated = hallucinated

    def predict_single(
        self, product_smiles: str, reactants_smiles: str
    ) -> dict[str, bool]:
        return {"is_hallucination": reactants_smiles in self.hallucinated}


class TestBuildMlChecker:
    def test_retains_clean_aspirin_precursors(self) -> None:
        classifier = StubClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID], [ACETIC_ACID]])

        assert status == 200
        assert retained == [[SALICYLIC_ACID], [ACETIC_ACID]]
        assert classifier.calls == [
            (ASPIRIN, SALICYLIC_ACID),
            (ASPIRIN, ACETIC_ACID),
        ]

    def test_removes_hallucinated_pathway(self) -> None:
        classifier = SelectiveClassifier(hallucinated={SALICYLIC_ACID})
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID], [ACETIC_ACID]])

        assert status == 200
        assert retained == [[ACETIC_ACID]]

    def test_skips_invalid_smiles_without_calling_classifier(self) -> None:
        classifier = StubClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [["not>>>valid"]])

        assert status == 200
        assert retained == []
        assert classifier.calls == []

    def test_joins_multi_reactant_pathway_with_dot(self) -> None:
        classifier = StubClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(PARACETAMOL, [[PARA_AMINOPHENOL, ACETIC_ACID]])

        assert status == 200
        assert retained == [[PARA_AMINOPHENOL, ACETIC_ACID]]
        assert classifier.calls == [
            (PARACETAMOL, f"{PARA_AMINOPHENOL}.{ACETIC_ACID}"),
        ]

    def test_missing_key_defaults_to_hallucinated(self) -> None:
        """Classifier returning no is_hallucination key is treated as hallucinated."""
        classifier = StubClassifier(responses=[{"confidence": 0.5}])
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID]])

        assert status == 200
        assert retained == []

    def test_none_value_treated_as_hallucinated(self) -> None:
        """is_hallucination=None is not the same as False."""
        classifier = StubClassifier(responses=[{"is_hallucination": None}])
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID]])

        assert status == 200
        assert retained == []

    def test_normalizes_string_pathway_to_list(self) -> None:
        """A pathway given as a dot-joined string is normalized to a list."""
        classifier = StubClassifier()
        checker = build_ml_checker(classifier)
        combined = f"{SALICYLIC_ACID}.{ACETIC_ACID}"

        status, retained = checker(ASPIRIN, [combined])

        assert status == 200
        assert retained == [[combined]]
        assert classifier.calls == [(ASPIRIN, combined)]

    def test_skips_invalid_pathway_but_keeps_valid_ones(self) -> None:
        """An invalid pathway is skipped without dropping the valid ones."""
        classifier = StubClassifier()
        checker = build_ml_checker(classifier)

        status, retained = checker(ASPIRIN, [["not>>>valid"], [SALICYLIC_ACID]])

        assert status == 200
        assert retained == [[SALICYLIC_ACID]]
        assert classifier.calls == [(ASPIRIN, SALICYLIC_ACID)]

    def test_rejects_classifier_without_predict_single(self) -> None:
        with pytest.raises(TypeError, match="predict_single"):
            build_ml_checker(object())  # type: ignore[arg-type]


class TestFilterWithChecker:
    def test_retains_aligned_metadata_after_filtering(self) -> None:
        def keep_second(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 200, [pathways[1]]

        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN,
            [[SALICYLIC_ACID], [ACETIC_ACID]],
            ["hydrolysis", "deacetylation"],
            [0.3, 0.9],
            keep_second,
        )

        assert pathways == [[ACETIC_ACID]]
        assert explanations == ["deacetylation"]
        assert confidence == [0.9]

    def test_none_checker_passes_through(self) -> None:
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
        def fail(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 500, pathways

        result = filter_with_checker(ASPIRIN, [[SALICYLIC_ACID]], ["only"], [0.8], fail)

        assert result == ([], [], [])

    def test_duplicate_pathways_consume_metadata_in_order(self) -> None:
        def return_dupes(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 200, [[SALICYLIC_ACID], [SALICYLIC_ACID]]

        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN,
            [[SALICYLIC_ACID], [SALICYLIC_ACID], [ACETIC_ACID]],
            ["first", "second", "third"],
            [0.1, 0.2, 0.3],
            return_dupes,
        )

        assert pathways == [[SALICYLIC_ACID], [SALICYLIC_ACID]]
        assert explanations == ["first", "second"]
        assert confidence == [0.1, 0.2]

    def test_empty_pathways_returns_empty(self) -> None:
        pathways, explanations, confidence = filter_with_checker(
            ASPIRIN, [], [], [], None
        )
        assert pathways == []
        assert explanations == []
        assert confidence == []

    def test_unknown_pathway_from_checker_raises(self) -> None:
        """Checker returning a pathway absent from the original input is an error."""

        def inject_unknown(
            product: str, pathways: list[list[str]]
        ) -> tuple[int, list[list[str]]]:
            return 200, [[PARACETAMOL]]

        with pytest.raises(ValueError, match="not present in original input"):
            filter_with_checker(
                ASPIRIN, [[SALICYLIC_ACID]], ["hydrolysis"], [0.9], inject_unknown
            )

    def test_mismatched_lengths_raise(self) -> None:
        """Pathways and metadata must have the same length."""
        with pytest.raises(ValueError, match="same length"):
            filter_with_checker(
                ASPIRIN,
                [[SALICYLIC_ACID], [ACETIC_ACID]],
                ["only one"],
                [0.1, 0.2],
                None,
            )

    def test_shorter_confidence_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            filter_with_checker(
                ASPIRIN,
                [[SALICYLIC_ACID]],
                ["explanation"],
                [],
                None,
            )


class TestResolveHallucination:
    def test_none_mode(self) -> None:
        assert resolve_hallucination("none", None) is None

    def test_heuristic_mode(self) -> None:
        assert resolve_hallucination("heuristic", None) is None

    def test_case_insensitive(self) -> None:
        assert resolve_hallucination("NONE", None) is None
        assert resolve_hallucination("Heuristic", None) is None

    def test_ml_mode_without_classifier_raises(self) -> None:
        with pytest.raises(ValueError, match="requires"):
            resolve_hallucination("ml", None)

    def test_unknown_mode_raises(self) -> None:
        with pytest.raises(ValueError, match="hallucination_mode"):
            resolve_hallucination("unknown", None)

    def test_ml_mode_returns_working_checker(self) -> None:
        """ML mode builds a checker that actually filters hallucinated pathways."""
        classifier = SelectiveClassifier(hallucinated={SALICYLIC_ACID})
        checker = resolve_hallucination("ml", classifier)
        assert checker is not None

        status, retained = checker(ASPIRIN, [[SALICYLIC_ACID], [ACETIC_ACID]])

        assert status == 200
        assert retained == [[ACETIC_ACID]]
