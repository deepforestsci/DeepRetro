"""Human-auditable JSON report and provenance-sidecar generation."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from deepretro.evaluation.aggregate import BenchmarkEvaluation, ScoreVariant
from deepretro.evaluation.io import file_sha256
from deepretro.evaluation.types import ChemCensorScore

REPORT_SCHEMA_VERSION = 1
SCORE_PROVENANCE_SCHEMA_VERSION = 1


def write_benchmark_report(
    path: str | Path,
    evaluation: BenchmarkEvaluation,
    *,
    score_variant: ScoreVariant,
    predictions_path: str | Path,
    targets_path: str | Path,
    scores_path: str | Path,
) -> None:
    """Write aggregate metrics, diagnostics, and artifact hashes as JSON."""
    payload = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "score_variant": score_variant,
        "metrics": {
            "target_count": evaluation.target_count,
            "prediction_coverage": evaluation.prediction_coverage,
            "average_pt_max_cc": evaluation.average_pt_max_cc,
            "average_pt_top_k_cc": evaluation.average_pt_top_k_cc,
            "reference_target_count": evaluation.reference_target_count,
            "average_exact_match_at_k": evaluation.average_exact_match_at_k,
        },
        "targets": [asdict(target) for target in evaluation.target_evaluations],
        "artifacts": {
            "predictions_sha256": file_sha256(predictions_path),
            "targets_sha256": file_sha256(targets_path),
            "scores_sha256": file_sha256(scores_path),
        },
        "limitations": [
            "ChemCensor is precedent-based: a no-precedent score is not proof "
            "that a reaction is chemically invalid.",
            "PT metrics evaluate single-step proposals; they are not a complete "
            "multi-step route feasibility score.",
        ],
    }
    _write_json(path, payload)


def write_score_provenance(
    path: str | Path,
    scores: Sequence[ChemCensorScore],
    *,
    predictions_path: str | Path,
    targets_path: str | Path,
) -> None:
    """Write scorer/database identity and source-input hashes beside score JSONL."""
    status_counts = Counter(score.status for score in scores)
    provenance_rows = {
        _provenance_key(score): asdict(score.provenance)
        for score in scores
        if score.provenance is not None
    }
    payload = {
        "schema_version": SCORE_PROVENANCE_SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "artifacts": {
            "predictions_sha256": file_sha256(predictions_path),
            "targets_sha256": file_sha256(targets_path),
        },
        "score_status_counts": dict(sorted(status_counts.items())),
        "scorer_provenance": list(provenance_rows.values()),
    }
    _write_json(path, payload)


def _provenance_key(score: ChemCensorScore) -> tuple[object, ...]:
    """Return a hashable identity for deterministic provenance de-duplication."""
    assert score.provenance is not None
    provenance = score.provenance
    return (
        provenance.scorer_name,
        provenance.package_version,
        provenance.database_version,
        provenance.database_sha256,
        provenance.max_center_type,
        provenance.find_exact_match,
    )


def _write_json(path: str | Path, payload: object) -> None:
    """Write formatted JSON atomically enough for local offline artifacts."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
