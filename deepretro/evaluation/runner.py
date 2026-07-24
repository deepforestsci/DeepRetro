"""Orchestration for optional direct, local official scoring runs."""

from __future__ import annotations

from dataclasses import replace
from typing import Protocol, Sequence

from deepretro.evaluation.canonical import canonicalize_prediction
from deepretro.evaluation.types import (
    CanonicalPrediction,
    ChemCensorScore,
    SingleStepPrediction,
)


class PredictionScorer(Protocol):
    """Minimal interface shared by the official scorer and test doubles."""

    def score_prediction(self, prediction: CanonicalPrediction) -> ChemCensorScore:
        """Return a normalized score for one canonicalized prediction."""


def score_predictions(
    predictions: Sequence[SingleStepPrediction],
    scorer: PredictionScorer,
) -> list[ChemCensorScore]:
    """Score predictions once per canonical reaction in stable input order.

    The output deliberately remains one row per original prediction so that an
    aggregate can validate the score artifact against every input prediction. A
    duplicate's score is rebound to its own ID without repeating the expensive
    external precedent lookup.

    Prediction-id uniqueness is enforced at the I/O boundary
    (:func:`deepretro.evaluation.io.load_predictions` on input and
    :func:`deepretro.evaluation.io.write_scores` on output), so it is not
    re-checked here.
    """
    scores: list[ChemCensorScore] = []
    by_reaction: dict[str, ChemCensorScore] = {}
    for prediction in predictions:
        canonical_prediction = canonicalize_prediction(prediction)
        reaction_smiles = canonical_prediction.canonical_reaction_smiles
        if reaction_smiles is not None and reaction_smiles in by_reaction:
            scores.append(
                replace(
                    by_reaction[reaction_smiles],
                    prediction_id=prediction.prediction_id,
                )
            )
            continue

        score = scorer.score_prediction(canonical_prediction)
        _validate_scorer_output(canonical_prediction, score)
        scores.append(score)
        if reaction_smiles is not None:
            by_reaction[reaction_smiles] = score
    return scores


def _validate_scorer_output(
    prediction: CanonicalPrediction,
    score: ChemCensorScore,
) -> None:
    """Prevent a third-party scorer bug from producing a mismatched artifact."""
    prediction_id = prediction.prediction.prediction_id
    if score.prediction_id != prediction_id:
        raise ValueError(
            f"scorer returned prediction_id {score.prediction_id!r} for {prediction_id!r}"
        )
    if prediction.is_valid:
        if score.canonical_reaction_smiles != prediction.canonical_reaction_smiles:
            raise ValueError(
                f"scorer returned a reaction different from prediction {prediction_id!r}"
            )
    elif score.status != "invalid":
        raise ValueError(
            f"scorer returned non-invalid status for invalid prediction {prediction_id!r}"
        )
