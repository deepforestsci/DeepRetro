"""Bootstrap a small subtagged ground-truth seed from the merged dataset."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from deepretro.data.subtags import (
    PRIMARY_TO_SECONDARY,
    bootstrap_seed_rows,
    summarize_label_counts_by_field,
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_report(rows: list[dict[str, str]]) -> dict[str, object]:
    return {
        "row_count": len(rows),
        "label_counts_by_source": summarize_label_counts_by_field(rows, "source"),
        "label_counts_by_primary_tag": summarize_label_counts_by_field(
            rows, "primary_hallucination_tag"
        ),
        "label_counts_by_secondary_tag": summarize_label_counts_by_field(
            rows, "secondary_hallucination_tag"
        ),
        "taxonomy": {
            primary: list(secondary)
            for primary, secondary in PRIMARY_TO_SECONDARY.items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default=Path("data/hallucination_dataset_merged.csv"),
        type=Path,
    )
    parser.add_argument(
        "--output",
        default=Path("data/hallucination_ground_truth_seed.csv"),
        type=Path,
    )
    parser.add_argument(
        "--report",
        default=Path("data/hallucination_ground_truth_seed_report.json"),
        type=Path,
    )
    parser.add_argument("--max-positive", type=int, default=20)
    parser.add_argument("--max-clean", type=int, default=10)
    args = parser.parse_args()

    rows = read_rows(args.input)
    seed_rows = bootstrap_seed_rows(
        rows,
        max_positive=args.max_positive,
        max_clean=args.max_clean,
    )
    write_rows(args.output, seed_rows)
    args.report.write_text(json.dumps(build_report(seed_rows), indent=2), encoding="utf-8")
    print(f"Wrote {len(seed_rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
