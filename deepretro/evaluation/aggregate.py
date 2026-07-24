"""Per-target and benchmark ChemCensor plausibility/diversity aggregates."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Literal, Sequence

from deepretro.evaluation.canonical import canonicalize_prediction, canonicalize_smiles
from deepretro.evaluation.types import (
    ChemCensorScore,
    EvaluationTarget,
    SingleStepPrediction,
)

ScoreVariant = Literal["with_functional_groups", "without_functional_groups"]


@dataclass(frozen=True)
class TargetEvaluation:
    """All reportable metrics for one benchmark target."""

    target_id: str
    k: int
    pt_max_cc: float
    pt_top_k_cc: float
    top_k_values: tuple[float, ...]
    exact_match_at_k: bool | None
    unique_prediction_count: int
    duplicate_prediction_count: int
    invalid_prediction_count: int
    missing_score_count: int
    processing_failed_count: int
    no_precedent_count: int
    scored_prediction_count: int


@dataclass(frozen=True)
class BenchmarkEvaluation:
    """Mean target-level metrics with denominator and coverage diagnostics."""

    target_evaluations: tuple[TargetEvaluation, ...]
    target_count: int
    prediction_coverage: float
    average_pt_max_cc: float
    average_pt_top_k_cc: float
    reference_target_count: int
    average_exact_match_at_k: float | None


def evaluate_target(
    target: EvaluationTarget,
    predictions: Sequence[SingleStepPrediction],
    scores: Sequence[ChemCensorScore],
    *,
    k: int,
    score_variant: ScoreVariant = "with_functional_groups",
) -> TargetEvaluation:
    """Compute PT-Max and zero-padded PT-Top-K for one explicit target.

    The maximum examines every unique valid prediction supplied for the target.
    PT-Top-K instead follows ranked output slots: invalid output consumes a zero
    slot, while a canonical duplicate consumes no slot and cannot manufacture
    diversity.  Missing/failed external scores are zero in aggregate values but
    have separate diagnostic counts.
    """
    _validate_metric_arguments(k, score_variant)
    target_product = canonicalize_smiles(target.product_smiles)
    if target_product is None:
        raise ValueError(f"target {target.target_id!r} has invalid product SMILES")

    _validate_prediction_records(target, predictions)
    score_by_id = _index_scores(scores, predictions)
    reference_reaction = _canonical_reference_reaction(target)

    unique_reactions: set[str] = set()
    top_k_values: list[float] = []
    top_k_reactions: list[str | None] = []
    pt_max_cc = 0.0
    duplicate_prediction_count = 0
    invalid_prediction_count = 0
    missing_score_count = 0
    processing_failed_count = 0
    no_precedent_count = 0
    scored_prediction_count = 0

    ranked_predictions = sorted(enumerate(predictions), key=lambda item: (item[1].rank, item[0]))
    for _input_index, prediction in ranked_predictions:
        canonical_prediction = canonicalize_prediction(prediction)
        score = score_by_id.get(prediction.prediction_id)
        if not canonical_prediction.is_valid:
            invalid_prediction_count += 1
            _validate_invalid_score(score, prediction.prediction_id)
            if len(top_k_values) < k:
                top_k_values.append(0.0)
                top_k_reactions.append(None)
            continue

        assert canonical_prediction.canonical_product_smiles is not None
        if canonical_prediction.canonical_product_smiles != target_product:
            raise ValueError(
                f"prediction {prediction.prediction_id!r} product does not match "
                f"target {target.target_id!r}"
            )
        assert canonical_prediction.canonical_reaction_smiles is not None
        _validate_valid_score(score, canonical_prediction.canonical_reaction_smiles)
        reaction_smiles = canonical_prediction.canonical_reaction_smiles
        if reaction_smiles in unique_reactions:
            duplicate_prediction_count += 1
            continue
        unique_reactions.add(reaction_smiles)

        value, score_state = _selected_score_value(score, score_variant)
        if score_state == "missing":
            missing_score_count += 1
        elif score_state == "processing_failed":
            processing_failed_count += 1
        elif score_state == "no_precedent":
            no_precedent_count += 1
        else:
            scored_prediction_count += 1
        pt_max_cc = max(pt_max_cc, value)
        if len(top_k_values) < k:
            top_k_values.append(value)
            top_k_reactions.append(reaction_smiles)

    padded_values = tuple(top_k_values + [0.0] * (k - len(top_k_values)))
    exact_match_at_k = (
        None
        if reference_reaction is None
        else reference_reaction in top_k_reactions
    )
    return TargetEvaluation(
        target_id=target.target_id,
        k=k,
        pt_max_cc=pt_max_cc,
        pt_top_k_cc=sum(padded_values) / k,
        top_k_values=padded_values,
        exact_match_at_k=exact_match_at_k,
        unique_prediction_count=len(unique_reactions),
        duplicate_prediction_count=duplicate_prediction_count,
        invalid_prediction_count=invalid_prediction_count,
        missing_score_count=missing_score_count,
        processing_failed_count=processing_failed_count,
        no_precedent_count=no_precedent_count,
        scored_prediction_count=scored_prediction_count,
    )


def aggregate_benchmark(
    targets: Sequence[EvaluationTarget],
    predictions: Sequence[SingleStepPrediction],
    scores: Sequence[ChemCensorScore],
    *,
    k: int,
    score_variant: ScoreVariant = "with_functional_groups",
) -> BenchmarkEvaluation:
    """Average target-level metrics over the full explicit target denominator."""
    _validate_metric_arguments(k, score_variant)
    targets_by_id: dict[str, EvaluationTarget] = {}
    for target in targets:
        if target.target_id in targets_by_id:
            raise ValueError(f"duplicate target_id: {target.target_id!r}")
        targets_by_id[target.target_id] = target
    if not targets_by_id:
        raise ValueError("targets must contain at least one target")

    predictions_by_target: dict[str, list[SingleStepPrediction]] = defaultdict(list)
    prediction_ids: set[str] = set()
    for prediction in predictions:
        if prediction.prediction_id in prediction_ids:
            raise ValueError(f"duplicate prediction_id: {prediction.prediction_id!r}")
        prediction_ids.add(prediction.prediction_id)
        if prediction.target_id not in targets_by_id:
            raise ValueError(
                f"prediction {prediction.prediction_id!r} references unknown target "
                f"{prediction.target_id!r}"
            )
        predictions_by_target[prediction.target_id].append(prediction)

    score_by_id: dict[str, ChemCensorScore] = {}
    for score in scores:
        if score.prediction_id not in prediction_ids:
            raise ValueError(
                f"score references unknown prediction_id: {score.prediction_id!r}"
            )
        if score.prediction_id in score_by_id:
            raise ValueError(f"duplicate score prediction_id: {score.prediction_id!r}")
        score_by_id[score.prediction_id] = score

    per_target = tuple(
        evaluate_target(
            target,
            predictions_by_target.get(target.target_id, []),
            [
                score_by_id[prediction.prediction_id]
                for prediction in predictions_by_target.get(target.target_id, [])
                if prediction.prediction_id in score_by_id
            ],
            k=k,
            score_variant=score_variant,
        )
        for target in targets
    )
    exact_matches = [
        result.exact_match_at_k
        for result in per_target
        if result.exact_match_at_k is not None
    ]
    target_count = len(per_target)
    return BenchmarkEvaluation(
        target_evaluations=per_target,
        target_count=target_count,
        prediction_coverage=(
            sum(bool(predictions_by_target[target.target_id]) for target in targets)
            / target_count
        ),
        average_pt_max_cc=sum(result.pt_max_cc for result in per_target) / target_count,
        average_pt_top_k_cc=(
            sum(result.pt_top_k_cc for result in per_target) / target_count
        ),
        reference_target_count=len(exact_matches),
        average_exact_match_at_k=(
            None
            if not exact_matches
            else sum(bool(value) for value in exact_matches) / len(exact_matches)
        ),
    )


def _validate_metric_arguments(k: int, score_variant: str) -> None:
    """Keep aggregate denominator and score selection unambiguous."""
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")
    if score_variant not in {"with_functional_groups", "without_functional_groups"}:
        raise ValueError(
            "score_variant must be 'with_functional_groups' or "
            "'without_functional_groups'"
        )


def _validate_prediction_records(
    target: EvaluationTarget,
    predictions: Sequence[SingleStepPrediction],
) -> None:
    """Reject ambiguous rank/ID inputs before scoring a target."""
    prediction_ids: set[str] = set()
    ranks: set[int] = set()
    for prediction in predictions:
        if prediction.target_id != target.target_id:
            raise ValueError(
                f"prediction {prediction.prediction_id!r} belongs to "
                f"{prediction.target_id!r}, not {target.target_id!r}"
            )
        if prediction.prediction_id in prediction_ids:
            raise ValueError(f"duplicate prediction_id: {prediction.prediction_id!r}")
        if prediction.rank in ranks:
            raise ValueError(
                f"target {target.target_id!r} has duplicate rank {prediction.rank}"
            )
        prediction_ids.add(prediction.prediction_id)
        ranks.add(prediction.rank)


def _index_scores(
    scores: Sequence[ChemCensorScore],
    predictions: Sequence[SingleStepPrediction],
) -> dict[str, ChemCensorScore]:
    """Index scores while rejecting stale and duplicate artifact rows."""
    prediction_ids = {prediction.prediction_id for prediction in predictions}
    indexed: dict[str, ChemCensorScore] = {}
    for score in scores:
        if score.prediction_id not in prediction_ids:
            raise ValueError(
                f"score references unknown prediction_id: {score.prediction_id!r}"
            )
        if score.prediction_id in indexed:
            raise ValueError(f"duplicate score prediction_id: {score.prediction_id!r}")
        indexed[score.prediction_id] = score
    return indexed


def _validate_invalid_score(score: ChemCensorScore | None, prediction_id: str) -> None:
    """Ensure an artifact cannot claim a score for an invalid current reaction."""
    if score is not None and score.status != "invalid":
        raise ValueError(
            f"invalid prediction {prediction_id!r} has a non-invalid score artifact"
        )


def _validate_valid_score(
    score: ChemCensorScore | None,
    canonical_reaction_smiles: str,
) -> None:
    """Reject score records produced for a different canonical prediction."""
    if score is None:
        return
    if score.canonical_reaction_smiles != canonical_reaction_smiles:
        raise ValueError(
            "score artifact reaction does not match current canonical prediction: "
            f"{score.prediction_id!r}"
        )


def _selected_score_value(
    score: ChemCensorScore | None,
    score_variant: ScoreVariant,
) -> tuple[float, str]:
    """Normalize unavailable/failed results to zero without concealing status."""
    if score is None:
        return 0.0, "missing"
    if score.status == "processing_failed":
        return 0.0, "processing_failed"
    if score.status == "invalid":
        return 0.0, "missing"
    if score.status == "no_precedent":
        return 0.0, "no_precedent"
    value = getattr(score, score_variant)
    assert value is not None
    return float(value), "scored"


def _canonical_reference_reaction(target: EvaluationTarget) -> str | None:
    """Canonicalize an optional literature precursor set using the same rules."""
    if target.reference_precursor_smiles is None:
        return None
    reference = SingleStepPrediction(
        prediction_id=f"{target.target_id}:reference",
        target_id=target.target_id,
        rank=1,
        product_smiles=target.product_smiles,
        precursor_smiles=target.reference_precursor_smiles,
    )
    canonical_reference = canonicalize_prediction(reference)
    if not canonical_reference.is_valid:
        raise ValueError(
            f"target {target.target_id!r} has invalid reference precursor SMILES"
        )
    return canonical_reference.canonical_reaction_smiles
