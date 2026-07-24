"""Versioned JSONL input/output for offline evaluation runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from deepretro.evaluation.types import (
    ChemCensorProvenance,
    ChemCensorScore,
    EvaluationTarget,
    SingleStepPrediction,
)

SCORE_ARTIFACT_SCHEMA_VERSION = 1
_PROVENANCE_FIELDS = {
    "scorer_name",
    "package_version",
    "database_filename",
    "database_version",
    "database_sha256",
    "database_size_bytes",
    "max_center_type",
    "find_exact_match",
}


def load_predictions(path: str | Path) -> list[SingleStepPrediction]:
    """Load ranked precursor-set predictions from the documented JSONL schema."""
    predictions: list[SingleStepPrediction] = []
    prediction_ids: set[str] = set()
    ranks_by_target: dict[str, set[int]] = {}
    for line_number, record in _read_jsonl(path):
        _require_exact_keys(
            record,
            required={
                "prediction_id",
                "target_id",
                "rank",
                "product_smiles",
                "precursor_smiles",
            },
            optional={"metadata"},
            source=f"{path}:{line_number}",
        )
        precursor_smiles = record["precursor_smiles"]
        if not isinstance(precursor_smiles, list):
            raise ValueError(f"{path}:{line_number}: precursor_smiles must be a list")
        prediction = SingleStepPrediction(
            prediction_id=record["prediction_id"],
            target_id=record["target_id"],
            rank=record["rank"],
            product_smiles=record["product_smiles"],
            precursor_smiles=tuple(precursor_smiles),
            metadata=_metadata_or_empty(record.get("metadata"), path, line_number),
        )
        if prediction.prediction_id in prediction_ids:
            raise ValueError(f"{path}:{line_number}: duplicate prediction_id")
        target_ranks = ranks_by_target.setdefault(prediction.target_id, set())
        if prediction.rank in target_ranks:
            raise ValueError(
                f"{path}:{line_number}: duplicate rank {prediction.rank} for "
                f"target {prediction.target_id!r}"
            )
        prediction_ids.add(prediction.prediction_id)
        target_ranks.add(prediction.rank)
        predictions.append(prediction)
    return predictions


def load_targets(path: str | Path) -> list[EvaluationTarget]:
    """Load the explicit benchmark target denominator from JSONL."""
    targets: list[EvaluationTarget] = []
    target_ids: set[str] = set()
    for line_number, record in _read_jsonl(path):
        _require_exact_keys(
            record,
            required={"target_id", "product_smiles"},
            optional={"reference_precursor_smiles", "metadata"},
            source=f"{path}:{line_number}",
        )
        reference_precursors = record.get("reference_precursor_smiles")
        if reference_precursors is not None and not isinstance(reference_precursors, list):
            raise ValueError(
                f"{path}:{line_number}: reference_precursor_smiles must be a list"
            )
        target = EvaluationTarget(
            target_id=record["target_id"],
            product_smiles=record["product_smiles"],
            reference_precursor_smiles=(
                None if reference_precursors is None else tuple(reference_precursors)
            ),
            metadata=_metadata_or_empty(record.get("metadata"), path, line_number),
        )
        if target.target_id in target_ids:
            raise ValueError(f"{path}:{line_number}: duplicate target_id")
        target_ids.add(target.target_id)
        targets.append(target)
    return targets


def write_scores(path: str | Path, scores: Sequence[ChemCensorScore]) -> None:
    """Write a deterministic, self-describing score artifact JSONL file."""
    score_ids: set[str] = set()
    for score in scores:
        if score.prediction_id in score_ids:
            raise ValueError(f"duplicate prediction_id: {score.prediction_id!r}")
        score_ids.add(score.prediction_id)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for score in scores:
            handle.write(json.dumps(_score_record(score), sort_keys=True))
            handle.write("\n")


def load_scores(path: str | Path) -> list[ChemCensorScore]:
    """Load a score artifact and validate its schema/status invariants."""
    scores: list[ChemCensorScore] = []
    score_ids: set[str] = set()
    for line_number, record in _read_jsonl(path):
        _require_exact_keys(
            record,
            required={
                "schema_version",
                "prediction_id",
                "canonical_reaction_smiles",
                "with_functional_groups",
                "without_functional_groups",
                "status",
                "provenance",
                "error",
            },
            optional=set(),
            source=f"{path}:{line_number}",
        )
        schema_version = record["schema_version"]
        if (
            isinstance(schema_version, bool)
            or not isinstance(schema_version, int)
            or schema_version != SCORE_ARTIFACT_SCHEMA_VERSION
        ):
            raise ValueError(
                f"{path}:{line_number}: schema_version must be integer "
                f"{SCORE_ARTIFACT_SCHEMA_VERSION}"
            )
        provenance_record = record["provenance"]
        if provenance_record is not None and not isinstance(provenance_record, Mapping):
            raise ValueError(f"{path}:{line_number}: provenance must be an object or null")
        if provenance_record is not None:
            _require_exact_keys(
                provenance_record,
                required=_PROVENANCE_FIELDS,
                optional=set(),
                source=f"{path}:{line_number}:provenance",
            )
        try:
            provenance = (
                None
                if provenance_record is None
                else ChemCensorProvenance(**dict(provenance_record))
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number}: invalid provenance: {exc}") from exc
        try:
            score = ChemCensorScore(
                prediction_id=record["prediction_id"],
                canonical_reaction_smiles=record["canonical_reaction_smiles"],
                with_functional_groups=record["with_functional_groups"],
                without_functional_groups=record["without_functional_groups"],
                status=record["status"],
                provenance=provenance,
                error=record["error"],
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number}: invalid score record: {exc}") from exc
        if score.prediction_id in score_ids:
            raise ValueError(f"{path}:{line_number}: duplicate prediction_id")
        score_ids.add(score.prediction_id)
        scores.append(score)
    return scores


def file_sha256(path: str | Path) -> str:
    """Return a streamed SHA-256 digest for artifact provenance sidecars."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_jsonl(path: str | Path) -> Iterable[tuple[int, dict[str, Any]]]:
    """Yield nonblank object records with file/line-contextual validation."""
    source = Path(path)
    with source.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            try:
                record = json.loads(raw_line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{source}:{line_number}: invalid JSON") from exc
            if not isinstance(record, dict):
                raise ValueError(f"{source}:{line_number}: each JSONL row must be an object")
            yield line_number, record


def _require_exact_keys(
    record: Mapping[str, Any],
    *,
    required: set[str],
    optional: set[str],
    source: str,
) -> None:
    """Reject missing and unknown fields instead of silently changing contracts."""
    missing = required - record.keys()
    unknown = record.keys() - required - optional
    if missing:
        raise ValueError(f"{source}: missing required fields: {sorted(missing)}")
    if unknown:
        raise ValueError(f"{source}: unknown fields: {sorted(unknown)}")


def _metadata_or_empty(value: Any, path: str | Path, line_number: int) -> Mapping[str, Any]:
    """Normalize optional metadata while keeping the primary schema strict."""
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"{path}:{line_number}: metadata must be an object")
    return dict(value)


def _score_record(score: ChemCensorScore) -> dict[str, Any]:
    """Serialize a score without relying on non-JSON dataclass internals."""
    return {
        "schema_version": SCORE_ARTIFACT_SCHEMA_VERSION,
        "prediction_id": score.prediction_id,
        "canonical_reaction_smiles": score.canonical_reaction_smiles,
        "with_functional_groups": score.with_functional_groups,
        "without_functional_groups": score.without_functional_groups,
        "status": score.status,
        "provenance": (
            None
            if score.provenance is None
            else {
                "scorer_name": score.provenance.scorer_name,
                "package_version": score.provenance.package_version,
                "database_filename": score.provenance.database_filename,
                "database_version": score.provenance.database_version,
                "database_sha256": score.provenance.database_sha256,
                "database_size_bytes": score.provenance.database_size_bytes,
                "max_center_type": score.provenance.max_center_type,
                "find_exact_match": score.provenance.find_exact_match,
            }
        ),
        "error": score.error,
    }
