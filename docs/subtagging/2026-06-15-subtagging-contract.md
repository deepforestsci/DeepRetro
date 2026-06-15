# Subtagging Contract

Date: 2026-06-15

## Purpose

This PR adds deterministic subtype enrichment for hallucination annotation rows.

## Canonical Module

- `deepretro/data/subtagging.py`

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
