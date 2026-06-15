# Subtagging Contract

Date: 2026-06-15

## Purpose

This PR adds deterministic subtype enrichment for hallucination annotation rows and a lightweight supervised train/eval/predict pipeline for subtype tags.

The dataset base for this work was split out from PR `#219` by taking only the merged hallucination dataset artifact.
The reviewed Ohuamine seed CSVs were then layered on top locally, while classifier and featurizer changes from `#219` were left out.

## Canonical Module

- `deepretro/data/subtagging.py`

## Pipeline Script

- `scripts/train_eval_predict_subtags.py`

The script supports `train`, `train-eval`, `eval`, and `predict` subcommands. It writes a JSON model artifact and appends `predicted_subtype_primary` during prediction.

## Canonical Output Columns

- `subtype_primary`
- `subtype_secondary`

## Input Expectations

The enrichment script expects CSV rows with:

- `label`
- either `category` or `reason`

Reviewed rows use `category` first when present. Free-text rows fall back to `reason`.

Rows without either field are still enriched:

- `label == 0` -> `subtype_primary=valid_clean`
- `label == 1` -> `subtype_primary=unclassified_hallucination`

## Out of Scope

This PR does not change model training, classifier thresholds, merged-dataset generation, or feature engineering.
