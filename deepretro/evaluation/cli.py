"""Command-line interface for optional scoring and offline aggregation."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from deepretro.evaluation.aggregate import aggregate_benchmark
from deepretro.evaluation.chemcensor import (
    ChemCensorUnavailableError,
    OfficialChemCensorScorer,
)
from deepretro.evaluation.io import load_predictions, load_scores, load_targets, write_scores
from deepretro.evaluation.report import write_benchmark_report, write_score_provenance
from deepretro.evaluation.runner import score_predictions


def main(argv: Sequence[str] | None = None) -> int:
    """Run ``score`` or ``aggregate`` and return a shell-compatible status."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "score":
            _run_score(args)
        elif args.command == "aggregate":
            _run_aggregate(args)
        else:  # pragma: no cover - argparse enforces the available choices
            parser.error(f"unknown command: {args.command}")
    except (ChemCensorUnavailableError, FileNotFoundError, ValueError) as exc:
        print(f"deepretro evaluation: {exc}", file=sys.stderr)
        return 2
    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Create the intentionally narrow offline-only command contract."""
    parser = argparse.ArgumentParser(
        prog="python -m deepretro.evaluation",
        description=(
            "Offline ChemCensor plausibility evaluation for DeepRetro. "
            "'aggregate' runs on any supported interpreter; 'score' additionally "
            "requires Python >=3.12 and a local ChemCensor install."
        ),
    )
    commands = parser.add_subparsers(dest="command", required=True)

    score = commands.add_parser(
        "score",
        help=(
            "run an approved local official ChemCensor installation "
            "(requires Python >=3.12 and a local ChemCensor install)"
        ),
    )
    _add_input_arguments(score)
    score.add_argument("--database", required=True, help="local official SQLite database")
    score.add_argument("--database-version", required=True)
    score.add_argument("--scores-out", required=True)
    score.add_argument("--provenance-out", required=True)
    score.add_argument("--cache", help="optional local SQLite score cache")
    score.add_argument("--max-center-type", type=int, default=4)
    score.add_argument(
        "--no-find-exact-match",
        action="store_true",
        help="disable official exact-reaction lookup optimization",
    )

    aggregate = commands.add_parser(
        "aggregate",
        help="aggregate an existing score artifact without ChemCensor installed",
    )
    _add_input_arguments(aggregate)
    aggregate.add_argument("--scores", required=True)
    aggregate.add_argument("--k", type=int, default=10)
    aggregate.add_argument(
        "--score-variant",
        choices=("with_functional_groups", "without_functional_groups"),
        default="with_functional_groups",
    )
    aggregate.add_argument("--report-out", required=True)
    return parser


def _add_input_arguments(parser: argparse.ArgumentParser) -> None:
    """Share the explicit benchmark contract between both subcommands."""
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--targets", required=True)


def _run_score(args: argparse.Namespace) -> None:
    """Score every prediction with a locally installed official engine."""
    predictions = load_predictions(args.predictions)
    targets = load_targets(args.targets)
    _validate_prediction_targets(predictions, {target.target_id for target in targets})
    scorer = OfficialChemCensorScorer(
        args.database,
        database_version=args.database_version,
        max_center_type=args.max_center_type,
        find_exact_match=not args.no_find_exact_match,
        cache_path=args.cache,
    )
    scores = score_predictions(predictions, scorer)
    write_scores(args.scores_out, scores)
    write_score_provenance(
        args.provenance_out,
        scores,
        predictions_path=args.predictions,
        targets_path=args.targets,
    )


def _run_aggregate(args: argparse.Namespace) -> None:
    """Compute benchmark metrics from a portable score JSONL artifact."""
    predictions = load_predictions(args.predictions)
    targets = load_targets(args.targets)
    scores = load_scores(args.scores)
    evaluation = aggregate_benchmark(
        targets,
        predictions,
        scores,
        k=args.k,
        score_variant=args.score_variant,
    )
    write_benchmark_report(
        args.report_out,
        evaluation,
        score_variant=args.score_variant,
        predictions_path=args.predictions,
        targets_path=args.targets,
        scores_path=args.scores,
    )


def _validate_prediction_targets(predictions, target_ids: set[str]) -> None:
    """Fail early when a score run would use an unknown target."""
    if not target_ids:
        raise ValueError("targets must contain at least one target")
    for prediction in predictions:
        if prediction.target_id not in target_ids:
            raise ValueError(
                f"prediction {prediction.prediction_id!r} references unknown target "
                f"{prediction.target_id!r}"
            )
