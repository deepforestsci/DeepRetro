"""Core-functionality tests for the offline ChemCensor evaluation pipeline.

These cover the main flow end to end -- canonicalization, score orchestration
with a scorer double, plausibility/diversity aggregation, and score-artifact
round-tripping -- rather than every helper and validation branch.
"""

from __future__ import annotations

from pathlib import Path

from deepretro.evaluation.aggregate import aggregate_benchmark, evaluate_target
from deepretro.evaluation.canonical import canonicalize_prediction
from deepretro.evaluation.io import load_scores, write_scores
from deepretro.evaluation.runner import score_predictions
from deepretro.evaluation.types import (
    CanonicalPrediction,
    ChemCensorProvenance,
    ChemCensorScore,
    EvaluationTarget,
    SingleStepPrediction,
)


def _provenance() -> ChemCensorProvenance:
    return ChemCensorProvenance(
        scorer_name="chemcensor",
        package_version="1.2.0",
        database_filename="ChemCensor.sqlite",
        database_version="U2-1.0.0",
        database_sha256="a" * 64,
        database_size_bytes=123,
        max_center_type=4,
        find_exact_match=True,
    )


class _FakeScorer:
    """Scorer double returning fixed dual scores and recording each call."""

    def __init__(self) -> None:
        self.calls: list[str] = []
        self.provenance = _provenance()

    def score_prediction(self, prediction: CanonicalPrediction) -> ChemCensorScore:
        self.calls.append(prediction.canonical_reaction_smiles or "invalid")
        if not prediction.is_valid:
            return ChemCensorScore(
                prediction_id=prediction.prediction.prediction_id,
                canonical_reaction_smiles=None,
                with_functional_groups=None,
                without_functional_groups=None,
                status="invalid",
                provenance=None,
                error="invalid",
            )
        return ChemCensorScore(
            prediction_id=prediction.prediction.prediction_id,
            canonical_reaction_smiles=prediction.canonical_reaction_smiles,
            with_functional_groups=3.0,
            without_functional_groups=4.0,
            status="scored",
            provenance=self.provenance,
        )


def test_canonicalization_builds_sorted_forward_reaction() -> None:
    """Precursors are canonicalized and order-normalized into precursors>>product."""
    prediction = SingleStepPrediction("p1", "target-1", 1, "OC(C)=O", ("O", "CCO"))

    canonical = canonicalize_prediction(prediction)

    assert canonical.is_valid
    assert canonical.canonical_reaction_smiles == "CCO.O>>CC(=O)O"


def test_canonicalization_flags_invalid_smiles() -> None:
    """Invalid model output is retained as an audit record, not discarded."""
    prediction = SingleStepPrediction("p1", "target-1", 1, "CC(=O)O", ("not-a-smiles",))

    canonical = canonicalize_prediction(prediction)

    assert not canonical.is_valid
    assert canonical.canonical_reaction_smiles is None
    assert canonical.error is not None


def test_score_predictions_scores_each_reaction_once_and_rebinds_ids() -> None:
    """Duplicate canonical reactions reuse one external lookup, one row per input."""
    scorer = _FakeScorer()
    predictions = [
        SingleStepPrediction("p1", "target-1", 1, "CC(=O)O", ("CCO", "O")),
        SingleStepPrediction("p2", "target-1", 2, "OC(C)=O", ("O", "C(O)C")),
    ]

    scores = score_predictions(predictions, scorer)

    assert scorer.calls == ["CCO.O>>CC(=O)O"]
    assert [score.prediction_id for score in scores] == ["p1", "p2"]
    assert [score.with_functional_groups for score in scores] == [3.0, 3.0]


def test_aggregate_rewards_distinct_predictions_over_full_target_denominator() -> None:
    """PT-Max sees all uniques; PT-Top-K pads slots; empty targets stay counted."""
    targets = [
        EvaluationTarget("target-1", "CC(=O)O", reference_precursor_smiles=("CCO", "O")),
        EvaluationTarget("target-2", "CCO"),
    ]
    predictions = [
        SingleStepPrediction("p1", "target-1", 1, "CC(=O)O", ("CCO", "O")),
        SingleStepPrediction("p2", "target-1", 2, "CC(=O)O", ("CCN", "O")),
        SingleStepPrediction("p3", "target-1", 3, "CC(=O)O", ("CCO", "O")),
    ]
    scores = [
        ChemCensorScore("p1", "CCO.O>>CC(=O)O", 3.0, 3.0, "scored", _provenance()),
        ChemCensorScore("p2", "CCN.O>>CC(=O)O", 4.0, 4.0, "scored", _provenance()),
        ChemCensorScore("p3", "CCO.O>>CC(=O)O", 3.0, 3.0, "scored", _provenance()),
    ]

    per_target = evaluate_target(targets[0], predictions, scores, k=3)
    assert per_target.pt_max_cc == 4.0
    assert per_target.unique_prediction_count == 2
    assert per_target.duplicate_prediction_count == 1
    assert per_target.exact_match_at_k is True

    benchmark = aggregate_benchmark(targets, predictions, scores, k=3)
    assert benchmark.target_count == 2
    assert benchmark.prediction_coverage == 0.5
    assert benchmark.average_pt_max_cc == 2.0  # (4.0 + 0.0) / 2 targets


def test_score_artifact_round_trips_through_jsonl(tmp_path: Path) -> None:
    """write_scores/load_scores preserve the dual score and its provenance."""
    scores = [
        ChemCensorScore("p1", "CCO.O>>CC(=O)O", 3.0, 4.0, "scored", _provenance()),
    ]
    artifact = tmp_path / "scores.jsonl"

    write_scores(artifact, scores)
    loaded = load_scores(artifact)

    assert loaded == scores
