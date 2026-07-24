"""End-to-end pipeline test for the ChemCensor evaluation subpackage.

Runs a small, self-contained benchmark through the full offline path -- score via
the public ``PredictionScorer`` protocol, write the JSONL artifact, then aggregate
through the real CLI entry point -- and asserts the resulting report. It exercises
every score state (scored, duplicate reaction, no_precedent, invalid) without the
optional ChemCensor engine, so it runs on any supported interpreter.

The ``score`` CLI subcommand is intentionally not exercised here: it requires the
external ChemCensor package (Python >=3.12 + approved database), which is not a test
dependency. Scoring is done in-process with a reference scorer instead.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from deepretro.evaluation.cli import main
from deepretro.evaluation.io import write_scores
from deepretro.evaluation.runner import score_predictions
from deepretro.evaluation.types import (
    CanonicalPrediction,
    ChemCensorProvenance,
    ChemCensorScore,
    SingleStepPrediction,
)

# A tiny benchmark exercising all four score states. `confidence` in metadata drives
# the reference scorer; 0.0 -> no_precedent, an unparseable precursor -> invalid, and
# a repeated canonical reaction -> duplicate.
_PREDICTIONS = [
    # target "acid": two distinct scored disconnections + one duplicate of the first.
    {"prediction_id": "a1", "target_id": "acid", "rank": 1,
     "product_smiles": "CC(=O)O", "precursor_smiles": ["CCO", "O"],
     "metadata": {"confidence": 0.9}},
    {"prediction_id": "a2", "target_id": "acid", "rank": 2,
     "product_smiles": "OC(C)=O", "precursor_smiles": ["CCO", "O"],  # dup reaction of a1
     "metadata": {"confidence": 0.5}},
    {"prediction_id": "a3", "target_id": "acid", "rank": 3,
     "product_smiles": "CC(=O)O", "precursor_smiles": ["CC(=O)Cl", "O"],
     "metadata": {"confidence": 0.7}},
    # target "ethanol": one no_precedent (0 confidence) + one invalid precursor.
    {"prediction_id": "b1", "target_id": "ethanol", "rank": 1,
     "product_smiles": "CCO", "precursor_smiles": ["CC=O", "[H][H]"],
     "metadata": {"confidence": 0.0}},
    {"prediction_id": "b2", "target_id": "ethanol", "rank": 2,
     "product_smiles": "CCO", "precursor_smiles": ["not-a-smiles"],
     "metadata": {"confidence": 0.5}},
]
_TARGETS = [
    {"target_id": "acid", "product_smiles": "CC(=O)O"},
    {"target_id": "ethanol", "product_smiles": "CCO"},
]


class _ConfidenceScorer:
    """Reference PredictionScorer mapping recorded confidence onto the score shape."""

    def __init__(self) -> None:
        self.provenance = ChemCensorProvenance(
            scorer_name="reference-confidence",
            package_version="0.1.0",
            database_filename="reference.none",
            database_version="ref-1",
            database_sha256=hashlib.sha256(b"reference").hexdigest(),
            database_size_bytes=9,
            max_center_type=4,
            find_exact_match=True,
        )

    def score_prediction(self, prediction: CanonicalPrediction) -> ChemCensorScore:
        pid = prediction.prediction.prediction_id
        if not prediction.is_valid:
            return ChemCensorScore(
                prediction_id=pid, canonical_reaction_smiles=None,
                with_functional_groups=None, without_functional_groups=None,
                status="invalid", provenance=None, error="invalid prediction",
            )
        confidence = float(prediction.prediction.metadata["confidence"])
        if confidence <= 0.0:
            value, status = 0.0, "no_precedent"
        else:
            value, status = round(1.0 + 4.0 * confidence, 3), "scored"
        return ChemCensorScore(
            prediction_id=pid,
            canonical_reaction_smiles=prediction.canonical_reaction_smiles,
            with_functional_groups=value, without_functional_groups=value,
            status=status, provenance=self.provenance,
        )


def _write_jsonl(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")


def test_full_pipeline_scores_and_aggregates_through_the_cli(tmp_path: Path) -> None:
    predictions_path = tmp_path / "predictions.jsonl"
    targets_path = tmp_path / "targets.jsonl"
    scores_path = tmp_path / "scores.jsonl"
    report_path = tmp_path / "report.json"
    _write_jsonl(predictions_path, _PREDICTIONS)
    _write_jsonl(targets_path, _TARGETS)

    # Score in-process via the public runner + protocol scorer, write the artifact.
    pred_objs = [
        SingleStepPrediction(
            prediction_id=p["prediction_id"], target_id=p["target_id"], rank=p["rank"],
            product_smiles=p["product_smiles"],
            precursor_smiles=tuple(p["precursor_smiles"]), metadata=p["metadata"],
        )
        for p in _PREDICTIONS
    ]
    scores = score_predictions(pred_objs, _ConfidenceScorer())
    statuses = sorted(s.status for s in scores)
    assert statuses == ["invalid", "no_precedent", "scored", "scored", "scored"]
    write_scores(scores_path, scores)

    # Aggregate through the real CLI entry point.
    exit_code = main([
        "aggregate",
        "--predictions", str(predictions_path),
        "--targets", str(targets_path),
        "--scores", str(scores_path),
        "--k", "3",
        "--report-out", str(report_path),
    ])
    assert exit_code == 0

    report = json.loads(report_path.read_text(encoding="utf-8"))
    metrics = report["metrics"]
    assert metrics["target_count"] == 2
    assert metrics["prediction_coverage"] == 1.0
    assert metrics["average_exact_match_at_k"] is None
    # acid: pt_max = 4.6; ethanol: pt_max = 0.0  -> mean 2.3
    assert metrics["average_pt_max_cc"] == pytest.approx(2.3)
    # acid pt_top3 = (4.6 + 3.8) / 3; ethanol = 0  -> mean of the two
    assert metrics["average_pt_top_k_cc"] == pytest.approx((8.4 / 3 + 0.0) / 2)

    targets = {t["target_id"]: t for t in report["targets"]}
    acid = targets["acid"]
    assert acid["pt_max_cc"] == pytest.approx(4.6)
    assert acid["top_k_values"] == [pytest.approx(4.6), pytest.approx(3.8), 0.0]
    assert acid["unique_prediction_count"] == 2
    assert acid["duplicate_prediction_count"] == 1
    assert acid["scored_prediction_count"] == 2
    assert acid["no_precedent_count"] == 0
    assert acid["invalid_prediction_count"] == 0

    ethanol = targets["ethanol"]
    assert ethanol["pt_max_cc"] == 0.0
    assert ethanol["no_precedent_count"] == 1
    assert ethanol["invalid_prediction_count"] == 1
    assert ethanol["scored_prediction_count"] == 0
