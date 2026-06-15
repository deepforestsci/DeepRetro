"""Add normalized subtype columns to a hallucination dataset CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from deepretro.data.subtagging import enrich_annotation_row


def enrich_dataset(input_path: Path, output_path: Path) -> None:
    """Read ``input_path`` and write an enriched CSV to ``output_path``."""

    with input_path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        fieldnames = list(reader.fieldnames or [])
        for name in ("subtype_primary", "subtype_secondary"):
            if name not in fieldnames:
                fieldnames.append(name)

        rows = [enrich_annotation_row(row) for row in reader]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Parse CLI arguments and enrich a CSV in place or to a new path."""

    parser = argparse.ArgumentParser(
        description="Add subtype_primary and subtype_secondary columns to a CSV.",
    )
    parser.add_argument("input_csv", type=Path, help="Path to the source CSV.")
    parser.add_argument("output_csv", type=Path, help="Path to the enriched CSV.")
    args = parser.parse_args()
    enrich_dataset(args.input_csv, args.output_csv)


if __name__ == "__main__":
    main()
