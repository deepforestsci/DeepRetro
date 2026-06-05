"""Compare dataset/model generations with repeated stratified splits.

This script answers four separate questions:

1. How did the old 15-feature model perform on the old dataset?
2. How does the current 39-feature model perform on the old dataset?
3. How does the same old 15-feature model perform on the new merged dataset?
4. How does the current 39-feature model perform on the new merged dataset?

Thresholds are selected on a validation slice inside each training fold, never
on the test fold.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem, RDLogger
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import StratifiedKFold

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from deepretro.data.subtags import annotate_subtags, summarize_label_counts_by_field


OLD15_INDICES = list(range(15))


COMPARISONS = [
    {
        "name": "old_data_old15",
        "dataset": "old",
        "feature_set": "old15",
        "description": "Old dataset with ECFP4 plus the original 15 domain features.",
    },
    {
        "name": "old_data_current39",
        "dataset": "old",
        "feature_set": "current39",
        "description": "Old dataset with ECFP6 plus the current 39 domain features.",
    },
    {
        "name": "new_data_old15",
        "dataset": "new",
        "feature_set": "old15",
        "description": "New merged dataset with the same old 15-feature model.",
    },
    {
        "name": "new_data_current39",
        "dataset": "new",
        "feature_set": "current39",
        "description": "New merged dataset with ECFP6 plus the current 39 domain features.",
    },
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    usable = []
    for row in rows:
        if row.get("label") not in {"0", "1"}:
            continue
        if Chem.MolFromSmiles(row.get("product", "")) is None:
            continue
        reactants = [part for part in row.get("reactants", "").split(".") if part]
        if not reactants or any(Chem.MolFromSmiles(part) is None for part in reactants):
            continue
        if row.get("source") == "expert" and (row.get("reason") or row.get("label") == "0"):
            row = annotate_subtags(row)
        usable.append(row)
    return usable


def dataset_audit(
    rows: list[dict[str, str]],
    slice_fields: list[str],
) -> dict[str, dict[str, dict[str, int]]]:
    return {
        field: summarize_label_counts_by_field(rows, field)
        for field in slice_fields
    }


def labels(rows: list[dict[str, str]]) -> np.ndarray:
    return np.array([int(row["label"]) for row in rows], dtype=int)


def old15_domain_matrix(rows: list[dict[str, str]]) -> np.ndarray:
    from deepretro.utils.domain_features import extract_domain_features_single

    return np.vstack(
        [
            extract_domain_features_single(row["product"], row["reactants"])[OLD15_INDICES]
            for row in rows
        ]
    )


def featurize(rows: list[dict[str, str]], feature_set: str) -> np.ndarray:
    from deepretro.featurizers.reactionstep import ReactionStepFeaturizer

    pairs = [(row["product"], row["reactants"]) for row in rows]
    if feature_set == "old15":
        fp = ReactionStepFeaturizer(
            radius=2,
            size=2048,
            use_domain_features=False,
        ).featurize(pairs)
        X = np.hstack([fp, old15_domain_matrix(rows)])
    elif feature_set == "current39":
        X = ReactionStepFeaturizer(
            radius=3,
            size=2048,
            use_domain_features=True,
        ).featurize(pairs)
    else:
        raise ValueError(f"Unknown feature set: {feature_set}")
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


def fit_model(X_train: np.ndarray, y_train: np.ndarray, seed: int) -> XGBClassifier:
    from xgboost import XGBClassifier

    scale_pos_weight = 1.0
    positives = int(np.sum(y_train == 1))
    negatives = int(np.sum(y_train == 0))
    if positives > 0:
        scale_pos_weight = negatives / positives

    model = XGBClassifier(
        max_depth=6,
        learning_rate=0.05,
        n_estimators=300,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_alpha=0.1,
        reg_lambda=2.0,
        min_child_weight=3,
        gamma=0.1,
        random_state=seed,
        eval_metric="auc",
        scale_pos_weight=scale_pos_weight,
    )
    model.fit(X_train, y_train, verbose=False)
    return model


def find_threshold(
    y_true: np.ndarray,
    probabilities: np.ndarray,
    metric: str,
) -> tuple[float, float]:
    thresholds = np.unique(np.r_[np.linspace(0.0, 1.0, 201), probabilities])
    best_threshold = 0.5
    best_score = -1.0
    for threshold in thresholds:
        pred = (probabilities >= threshold).astype(int)
        if metric == "accuracy":
            score = accuracy_score(y_true, pred)
        elif metric == "f1":
            score = f1_score(y_true, pred, zero_division=0)
        else:
            raise ValueError(f"Unknown threshold metric: {metric}")
        if score > best_score:
            best_threshold = float(threshold)
            best_score = float(score)
    return best_threshold, best_score


def metric_row(
    comparison: dict[str, str],
    fold: int,
    seed: int,
    y_test: np.ndarray,
    probabilities: np.ndarray,
    threshold: float,
    n_train: int,
    n_test: int,
) -> dict[str, Any]:
    pred = (probabilities >= threshold).astype(int)
    return {
        "comparison": comparison["name"],
        "dataset": comparison["dataset"],
        "feature_set": comparison["feature_set"],
        "fold": fold,
        "seed": seed,
        "n_train": n_train,
        "n_test": n_test,
        "test_positive_rate": float(np.mean(y_test)),
        "threshold": float(threshold),
        "roc_auc": float(roc_auc_score(y_test, probabilities)),
        "pr_auc": float(average_precision_score(y_test, probabilities)),
        "accuracy": float(accuracy_score(y_test, pred)),
        "precision": float(precision_score(y_test, pred, zero_division=0)),
        "recall": float(recall_score(y_test, pred, zero_division=0)),
        "f1": float(f1_score(y_test, pred, zero_division=0)),
    }


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    metrics = ["roc_auc", "pr_auc", "accuracy", "precision", "recall", "f1"]
    summary = []
    for comparison in [item["name"] for item in COMPARISONS]:
        group = [row for row in rows if row["comparison"] == comparison]
        if not group:
            continue
        out: dict[str, Any] = {
            "comparison": comparison,
            "dataset": group[0]["dataset"],
            "feature_set": group[0]["feature_set"],
            "folds": len(group),
            "n_train_mean": float(np.mean([row["n_train"] for row in group])),
            "n_test_mean": float(np.mean([row["n_test"] for row in group])),
            "test_positive_rate_mean": float(
                np.mean([row["test_positive_rate"] for row in group])
            ),
        }
        for metric in metrics:
            values = np.array([row[metric] for row in group], dtype=float)
            std = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            ci95 = float(1.96 * std / np.sqrt(len(values))) if len(values) > 1 else 0.0
            out[f"{metric}_mean"] = float(np.mean(values))
            out[f"{metric}_std"] = std
            out[f"{metric}_ci95"] = ci95
        summary.append(out)
    return summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, Any]:
    datasets = {
        "old": read_rows(args.old_data),
        "new": read_rows(args.new_data),
    }
    audits = {
        name: dataset_audit(rows, args.slice_fields)
        for name, rows in datasets.items()
    }
    dataset_labels = {name: labels(rows) for name, rows in datasets.items()}
    fold_rows = []

    if args.audit_only:
        report = {
            "old_data": str(args.old_data),
            "new_data": str(args.new_data),
            "old_rows_used": len(datasets["old"]),
            "new_rows_used": len(datasets["new"]),
            "dataset_audits": audits,
            "comparisons": COMPARISONS,
            "summary": [],
            "audit_only": True,
        }
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "generation_stats_report.json").write_text(
            json.dumps(report, indent=2),
            encoding="utf-8",
        )
        print(json.dumps(audits, indent=2))
        return report

    for comparison in COMPARISONS:
        rows = datasets[comparison["dataset"]]
        y = dataset_labels[comparison["dataset"]]
        if args.split_mode == "kfold":
            split = StratifiedKFold(
                n_splits=args.folds,
                shuffle=True,
                random_state=args.seed,
            )
            split_iter = split.split(np.zeros(len(y)), y)
        else:
            split = StratifiedShuffleSplit(
                n_splits=args.folds,
                train_size=args.train_size,
                random_state=args.seed,
            )
            split_iter = split.split(np.zeros(len(y)), y)

        for fold, (train_idx, test_idx) in enumerate(split_iter, start=1):
            train_rows = [rows[index] for index in train_idx]
            test_rows = [rows[index] for index in test_idx]
            y_train = labels(train_rows)
            y_test = labels(test_rows)

            inner_split = StratifiedShuffleSplit(
                n_splits=1,
                test_size=args.valid_size,
                random_state=args.seed + fold,
            )
            fit_local, valid_local = next(
                inner_split.split(np.zeros(len(y_train)), y_train)
            )

            X_train = featurize(train_rows, comparison["feature_set"])
            X_test = featurize(test_rows, comparison["feature_set"])
            model = fit_model(X_train[fit_local], y_train[fit_local], args.seed + fold)

            valid_prob = model.predict_proba(X_train[valid_local])[:, 1]
            threshold, _ = find_threshold(
                y_train[valid_local],
                valid_prob,
                args.threshold_metric,
            )
            test_prob = model.predict_proba(X_test)[:, 1]

            fold_rows.append(
                metric_row(
                    comparison=comparison,
                    fold=fold,
                    seed=args.seed + fold - 1,
                    y_test=y_test,
                    probabilities=test_prob,
                    threshold=threshold,
                    n_train=len(train_rows),
                    n_test=len(test_rows),
                )
            )

    summary = summarize(fold_rows)
    report = {
        "old_data": str(args.old_data),
        "new_data": str(args.new_data),
        "old_rows_used": len(datasets["old"]),
        "new_rows_used": len(datasets["new"]),
        "folds": args.folds,
        "split_mode": args.split_mode,
        "train_size": args.train_size,
        "valid_size_within_train": args.valid_size,
        "threshold_policy": f"{args.threshold_metric} threshold selected on an inner validation slice of each training fold.",
        "dataset_audits": audits,
        "comparisons": COMPARISONS,
        "summary": summary,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "fold_metrics.csv", fold_rows)
    (args.output_dir / "generation_stats_report.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-data", default=Path("data/hallucination_dataset.csv"), type=Path)
    parser.add_argument("--new-data", default=Path("data/hallucination_dataset_merged.csv"), type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument("--split-mode", choices=["shuffle", "kfold"], default="kfold")
    parser.add_argument("--train-size", type=float, default=0.70)
    parser.add_argument("--valid-size", type=float, default=0.20)
    parser.add_argument("--threshold-metric", choices=["f1", "accuracy"], default="f1")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument(
        "--slice-fields",
        nargs="*",
        default=["source", "target", "primary_hallucination_tag", "secondary_hallucination_tag"],
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
