# Subtagging Dataset Foundation

Date: 2026-06-15

## Purpose

This PR establishes the dataset inputs that the follow-up subtagging PR will operate on.

## Canonical Files

- `data/hallucination_dataset_merged.csv`
- `data/ohuamine_synthetic_hallucination_rows_20260611.csv`
- `data/ohuamine_synthetic_hallucination_rows_v3_20260611.csv`
- `docs/ohuamine_hallucination_categories_20260611.md`
- `docs/tetracyclic_review_outcome_20260611.md`

## Roles

- `hallucination_dataset_merged.csv`: latest merged base dataset from PR #219 lineage
- `ohuamine_synthetic_hallucination_rows_20260611.csv`: reviewed Ohuamine rows plus synthetic counterexamples
- `ohuamine_synthetic_hallucination_rows_v3_20260611.csv`: tighter Ohuamine v3 augmentation set
- `ohuamine_hallucination_categories_20260611.md`: human-written rationale for the reviewed Ohuamine categories
- `tetracyclic_review_outcome_20260611.md`: human-written rationale that tetracyclic review rows should not remain confirmed hallucinations

## Out of Scope

This PR does not change classifier code, featurizer defaults, domain features, or subtagging logic.
