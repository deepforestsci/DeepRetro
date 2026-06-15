"""Train, evaluate, and predict hallucination subtype tags from CSV rows."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

_SUBTAGGING_PATH = PROJECT_ROOT / "deepretro" / "data" / "subtagging.py"
_SUBTAGGING_SPEC = importlib.util.spec_from_file_location(
    "_deepretro_subtagging",
    _SUBTAGGING_PATH,
)
if _SUBTAGGING_SPEC is None or _SUBTAGGING_SPEC.loader is None:
    raise ImportError(f"Could not load subtagging module from {_SUBTAGGING_PATH}")
_SUBTAGGING_MODULE = importlib.util.module_from_spec(_SUBTAGGING_SPEC)
sys.modules[_SUBTAGGING_SPEC.name] = _SUBTAGGING_MODULE
_SUBTAGGING_SPEC.loader.exec_module(_SUBTAGGING_MODULE)

SubtypeTextClassifier = _SUBTAGGING_MODULE.SubtypeTextClassifier
train_eval_split = _SUBTAGGING_MODULE.train_eval_split


def load_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows as dictionaries."""

    with path.open("r", newline="", encoding="utf-8") as infile:
        return list(csv.DictReader(infile))


def train_model_file(input_csv: Path, model_json: Path) -> SubtypeTextClassifier:
    """Train a subtype classifier from ``input_csv`` and persist it."""

    classifier = SubtypeTextClassifier().fit(load_rows(input_csv))
    classifier.save(model_json)
    return classifier


def evaluate_model_file(input_csv: Path, model_json: Path) -> dict[str, object]:
    """Evaluate a persisted subtype classifier against ``input_csv``."""

    classifier = SubtypeTextClassifier.load(model_json)
    return classifier.evaluate(load_rows(input_csv))


def predict_file(input_csv: Path, model_json: Path, output_csv: Path) -> None:
    """Write ``input_csv`` with an added ``predicted_subtype_primary`` column."""

    rows = load_rows(input_csv)
    classifier = SubtypeTextClassifier.load(model_json)
    fieldnames = list(rows[0].keys()) if rows else []
    if "predicted_subtype_primary" not in fieldnames:
        fieldnames.append("predicted_subtype_primary")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row, prediction in zip(rows, classifier.predict(rows)):
            output_row = dict(row)
            output_row["predicted_subtype_primary"] = prediction
            writer.writerow(output_row)


def train_eval_file(
    input_csv: Path,
    model_json: Path,
    eval_fraction: float,
    seed: int,
) -> dict[str, object]:
    """Train on a deterministic split and evaluate on the held-out rows."""

    train_rows, eval_rows = train_eval_split(
        load_rows(input_csv),
        eval_fraction=eval_fraction,
        seed=seed,
    )
    classifier = SubtypeTextClassifier().fit(train_rows)
    classifier.save(model_json)
    return classifier.evaluate(eval_rows)


def main() -> None:
    """Parse CLI arguments for train/eval/predict subcommands."""

    parser = argparse.ArgumentParser(
        description="Train, evaluate, and predict hallucination subtype tags.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("input_csv", type=Path)
    train_parser.add_argument("model_json", type=Path)

    eval_parser = subparsers.add_parser("eval")
    eval_parser.add_argument("input_csv", type=Path)
    eval_parser.add_argument("model_json", type=Path)

    split_parser = subparsers.add_parser("train-eval")
    split_parser.add_argument("input_csv", type=Path)
    split_parser.add_argument("model_json", type=Path)
    split_parser.add_argument("--eval-fraction", type=float, default=0.2)
    split_parser.add_argument("--seed", type=int, default=42)

    predict_parser = subparsers.add_parser("predict")
    predict_parser.add_argument("input_csv", type=Path)
    predict_parser.add_argument("model_json", type=Path)
    predict_parser.add_argument("output_csv", type=Path)

    args = parser.parse_args()
    if args.command == "train":
        train_model_file(args.input_csv, args.model_json)
    elif args.command == "eval":
        print(json.dumps(evaluate_model_file(args.input_csv, args.model_json), indent=2))
    elif args.command == "train-eval":
        metrics = train_eval_file(
            args.input_csv,
            args.model_json,
            eval_fraction=args.eval_fraction,
            seed=args.seed,
        )
        print(json.dumps(metrics, indent=2))
    elif args.command == "predict":
        predict_file(args.input_csv, args.model_json, args.output_csv)


if __name__ == "__main__":
    main()
