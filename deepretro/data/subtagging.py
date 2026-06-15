"""Subtype tagging utilities for hallucination-dataset annotations.

This module converts free-text ``reason`` annotations into a small,
normalized subtype taxonomy that can be used for downstream analysis
without changing the current model-training behavior.
"""

from __future__ import annotations

import json
import math
import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

SUBTYPE_VALID_CLEAN = "valid_clean"
SUBTYPE_PROTECTING_GROUP_ISSUE = "protecting_group_issue"
SUBTYPE_NO_REAL_DISCONNECTION = "no_real_disconnection"
SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION = "wrong_starting_material_or_direction"
SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR = "atom_count_or_formula_error"
SUBTYPE_RING_TOPOLOGY_ERROR = "ring_topology_error"
SUBTYPE_REGIOCHEMICAL_ERROR = "regiochemical_error"
SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR = "bond_or_fragment_connectivity_error"
SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR = "unstable_intermediate_or_precursor"
SUBTYPE_UNCLASSIFIED_HALLUCINATION = "unclassified_hallucination"

ALL_SUBTYPES = (
    SUBTYPE_VALID_CLEAN,
    SUBTYPE_PROTECTING_GROUP_ISSUE,
    SUBTYPE_NO_REAL_DISCONNECTION,
    SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION,
    SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
    SUBTYPE_RING_TOPOLOGY_ERROR,
    SUBTYPE_REGIOCHEMICAL_ERROR,
    SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
    SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR,
    SUBTYPE_UNCLASSIFIED_HALLUCINATION,
)

SUBTYPE_TEXT_FIELDS = ("product", "reactants", "target", "source")
_ALL_SUBTYPE_SET = frozenset(ALL_SUBTYPES)
_NULLISH_REASONS = frozenset(("", "n/a", "na", "none", "null", "unknown", "-"))

_CATEGORY_SUBTYPE_MAP: dict[str, tuple[str, str | None]] = {
    "unsupported_skeletal_edit": (
        SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
        None,
    ),
    "protecting_group_swap_fragment_loss": (
        SUBTYPE_PROTECTING_GROUP_ISSUE,
        SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
    ),
    "missing_methylene_source": (
        SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
        None,
    ),
    "valid_hydrogenation": (SUBTYPE_VALID_CLEAN, None),
    "valid_boc_deprotection": (SUBTYPE_VALID_CLEAN, None),
    "valid_debenzylation": (SUBTYPE_VALID_CLEAN, None),
    "valid_lactam_closure": (SUBTYPE_VALID_CLEAN, None),
    "valid_acid_chloride_formation": (SUBTYPE_VALID_CLEAN, None),
    "valid_ester_hydrolysis": (SUBTYPE_VALID_CLEAN, None),
    "valid_reductive_amination": (SUBTYPE_VALID_CLEAN, None),
}


@dataclass(frozen=True)
class SubtypeRule:
    """Map one or more reason-text markers to a normalized subtype."""

    subtype: str
    markers: tuple[str, ...]


_PRIMARY_RULES: tuple[SubtypeRule, ...] = (
    SubtypeRule(
        SUBTYPE_REGIOCHEMICAL_ERROR,
        ("position hallucination", "wrong position", "regio"),
    ),
    SubtypeRule(
        SUBTYPE_PROTECTING_GROUP_ISSUE,
        (
            "unnecessary protective group",
            "unfavourable protective group",
            "protective group",
            "protecting group",
            "methoxy",
        ),
    ),
    SubtypeRule(
        SUBTYPE_NO_REAL_DISCONNECTION,
        (
            "didn't break into smaller",
            "did not break into smaller",
            "no real progress",
            "no chemical reaction seen",
            "core unchanged",
            "pure hallucination",
        ),
    ),
    SubtypeRule(
        SUBTYPE_WRONG_STARTING_MATERIAL_OR_DIRECTION,
        ("wrong starting material",),
    ),
    SubtypeRule(
        SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
        (
            "atom count",
            "carbon missing",
            "carbon count increased",
            "sudden change in atom count",
        ),
    ),
    SubtypeRule(
        SUBTYPE_RING_TOPOLOGY_ERROR,
        (
            "membered ring",
            "ring opening",
            "ring formation",
            "ring formed",
            "ring missing",
            "ring reduced",
        ),
    ),
    SubtypeRule(
        SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
        (
            "attachments hallucinated",
            "bond between",
            "covalency hallucination",
            "connecting carbons",
            "splits into two fragments",
            "fragment breaks into two parts",
        ),
    ),
    SubtypeRule(
        SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR,
        ("unstable",),
    ),
    SubtypeRule(
        SUBTYPE_VALID_CLEAN,
        (
            "default-clean",
            "overall chemistry looks good",
            "proper chemistry",
            "no mistakes",
        ),
    ),
)


_SECONDARY_RULES: tuple[SubtypeRule, ...] = (
    SubtypeRule(SUBTYPE_UNSTABLE_INTERMEDIATE_OR_PRECURSOR, ("unstable",)),
    SubtypeRule(
        SUBTYPE_ATOM_COUNT_OR_FORMULA_ERROR,
        ("atom count", "carbon missing", "carbon count increased"),
    ),
    SubtypeRule(
        SUBTYPE_RING_TOPOLOGY_ERROR,
        (
            "membered ring",
            "ring opening",
            "ring formation",
            "ring missing",
            "ring reduced",
        ),
    ),
    SubtypeRule(
        SUBTYPE_BOND_OR_FRAGMENT_CONNECTIVITY_ERROR,
        ("attachments hallucinated", "bond between", "connecting carbons"),
    ),
)


def _normalize_text(reason: str | None) -> str:
    """Return a case-normalized reason string."""

    return (reason or "").strip().lower()


def _match_subtype(text: str, rules: Iterable[SubtypeRule]) -> str | None:
    """Return the first matching subtype from ``rules``."""

    for rule in rules:
        if any(marker in text for marker in rule.markers):
            return rule.subtype
    return None


def _match_secondary_subtype(
    text: str,
    rules: Iterable[SubtypeRule],
    primary: str,
) -> str | None:
    """Return the first matching subtype that differs from ``primary``."""

    for rule in rules:
        if rule.subtype == primary:
            continue
        if any(marker in text for marker in rule.markers):
            return rule.subtype
    return None


def _assign_from_category(
    label: str | int,
    category: str | None,
) -> tuple[str, str | None] | None:
    """Return a subtype assignment from a reviewed category, if available."""

    category_text = _normalize_text(category)
    if not category_text:
        return None

    assignment = _CATEGORY_SUBTYPE_MAP.get(category_text)
    if assignment is not None:
        return assignment

    return None


def assign_subtypes(
    label: str | int,
    reason: str | None,
) -> tuple[str, str | None]:
    """Assign normalized subtypes from a row label and free-text reason.

    Parameters
    ----------
    label : str or int
        Binary hallucination label. ``0`` is treated as clean/valid and
        ``1`` as hallucinated.
    reason : str or None
        Free-text annotation describing why the row was labeled.

    Returns
    -------
    tuple[str, str | None]
        ``(subtype_primary, subtype_secondary)``.
    """

    text = _normalize_text(reason)
    label_text = str(label).strip()

    if text:
        primary = _match_subtype(text, _PRIMARY_RULES)
        if primary is None:
            primary = (
                SUBTYPE_UNCLASSIFIED_HALLUCINATION
                if label_text == "1"
                else SUBTYPE_VALID_CLEAN
            )

        secondary = _match_secondary_subtype(text, _SECONDARY_RULES, primary)
        return primary, secondary

    if label_text == "1":
        return SUBTYPE_UNCLASSIFIED_HALLUCINATION, None
    return SUBTYPE_VALID_CLEAN, None


def enrich_annotation_row(row: dict[str, str]) -> dict[str, str]:
    """Return a copy of ``row`` with subtype columns added."""

    enriched = dict(row)
    category_assignment = _assign_from_category(
        label=row.get("label", ""),
        category=row.get("category", ""),
    )
    if category_assignment is not None:
        primary, secondary = category_assignment
    else:
        primary, secondary = assign_subtypes(
            label=row.get("label", ""),
            reason=row.get("reason", ""),
        )
    enriched["subtype_primary"] = primary
    enriched["subtype_secondary"] = secondary or ""
    return enriched


def train_eval_split(
    rows: list[dict[str, str]],
    eval_fraction: float = 0.2,
    seed: int = 42,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Return deterministic train/eval row splits for subtype modelling."""

    if not rows:
        return [], []
    if not 0 < eval_fraction < 1:
        raise ValueError("eval_fraction must be between 0 and 1")

    shuffled = list(rows)
    random.Random(seed).shuffle(shuffled)
    eval_size = max(1, round(len(shuffled) * eval_fraction))
    if eval_size >= len(shuffled):
        eval_size = len(shuffled) - 1
    return shuffled[eval_size:], shuffled[:eval_size]


def _row_target_subtype(row: dict[str, str]) -> str:
    """Resolve the supervised target subtype for a labelled row."""

    subtype = row.get("subtype_primary", "").strip()
    if subtype:
        if subtype not in _ALL_SUBTYPE_SET:
            raise ValueError(f"Unknown subtype_primary: {subtype}")
        return subtype

    category = row.get("category", "").strip()
    if category:
        category_assignment = _CATEGORY_SUBTYPE_MAP.get(_normalize_text(category))
        if category_assignment is None:
            raise ValueError(f"Unknown subtype category: {category}")
        return category_assignment[0]

    reason = row.get("reason", "").strip()
    reason_text = _normalize_text(reason)
    if reason_text and reason_text not in _NULLISH_REASONS:
        reason_target = _match_subtype(reason_text, _PRIMARY_RULES)
        if reason_target is not None:
            return reason_target

    raise ValueError(
        "Rows used for subtype training/evaluation require a supervised subtype target "
        "via subtype_primary, category, or reason."
    )


def _row_text(row: dict[str, str]) -> str:
    """Build non-label text used as features for subtype prediction."""

    return " ".join(str(row.get(field, "")) for field in SUBTYPE_TEXT_FIELDS)


def _tokenize_for_subtype_model(text: str) -> list[str]:
    """Return simple lexical and character-ngram features."""

    normalized = text.strip().lower()
    if not normalized:
        return []

    word_tokens = re.findall(r"[a-z0-9@+\-\[\]\(\)=#/\\]+", normalized)
    compact = re.sub(r"\s+", "", normalized)
    char_tokens = [
        f"char:{compact[index:index + size]}"
        for size in (2, 3, 4)
        for index in range(max(0, len(compact) - size + 1))
    ]
    return word_tokens + char_tokens


@dataclass
class SubtypeTextClassifier:
    """Small multinomial Naive Bayes classifier for subtype tags."""

    label_counts: dict[str, int] | None = None
    token_counts_by_label: dict[str, dict[str, int]] | None = None
    token_totals_by_label: dict[str, int] | None = None
    vocabulary: set[str] | None = None
    alpha: float = 1.0

    def fit(self, rows: list[dict[str, str]]) -> "SubtypeTextClassifier":
        """Train on labelled CSV rows and return ``self``."""

        label_counts: Counter[str] = Counter()
        token_counts: dict[str, Counter[str]] = defaultdict(Counter)
        token_totals: Counter[str] = Counter()
        vocabulary: set[str] = set()

        for row in rows:
            label = _row_target_subtype(row)
            tokens = _tokenize_for_subtype_model(_row_text(row))
            label_counts[label] += 1
            token_counts[label].update(tokens)
            token_totals[label] += len(tokens)
            vocabulary.update(tokens)

        if not label_counts:
            raise ValueError("Cannot train subtype classifier without rows")

        self.label_counts = dict(label_counts)
        self.token_counts_by_label = {
            label: dict(counts) for label, counts in token_counts.items()
        }
        self.token_totals_by_label = dict(token_totals)
        self.vocabulary = vocabulary
        return self

    def predict_row(self, row: dict[str, str]) -> str:
        """Predict ``subtype_primary`` for one row."""

        if not self.label_counts:
            raise ValueError("SubtypeTextClassifier must be fitted before prediction")

        tokens = _tokenize_for_subtype_model(_row_text(row))
        total_rows = sum(self.label_counts.values())
        labels = sorted(self.label_counts)
        vocabulary_size = max(1, len(self.vocabulary or set()))
        best_label = labels[0]
        best_score = -math.inf

        for label in labels:
            label_count = self.label_counts[label]
            token_counts = (self.token_counts_by_label or {}).get(label, {})
            token_total = (self.token_totals_by_label or {}).get(label, 0)
            score = math.log(label_count / total_rows)
            denominator = token_total + self.alpha * vocabulary_size
            for token in tokens:
                token_count = token_counts.get(token, 0)
                score += math.log((token_count + self.alpha) / denominator)
            if score > best_score:
                best_label = label
                best_score = score

        return best_label

    def predict(self, rows: list[dict[str, str]]) -> list[str]:
        """Predict ``subtype_primary`` for multiple rows."""

        return [self.predict_row(row) for row in rows]

    def evaluate(self, rows: list[dict[str, str]]) -> dict[str, object]:
        """Evaluate subtype predictions against resolved row targets."""

        if not rows:
            return {"accuracy": 0.0, "correct": 0, "total": 0, "labels": []}

        y_true = [_row_target_subtype(row) for row in rows]
        y_pred = self.predict(rows)
        correct = sum(1 for expected, actual in zip(y_true, y_pred) if expected == actual)
        return {
            "accuracy": correct / len(rows),
            "correct": correct,
            "total": len(rows),
            "labels": sorted(set(y_true)),
        }

    def save(self, path: str | Path) -> None:
        """Persist the trained classifier as JSON."""

        if not self.label_counts:
            raise ValueError("Cannot save an unfitted subtype classifier")

        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "label_counts": self.label_counts,
            "token_counts_by_label": self.token_counts_by_label or {},
            "token_totals_by_label": self.token_totals_by_label or {},
            "vocabulary": sorted(self.vocabulary or set()),
            "alpha": self.alpha,
        }
        output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "SubtypeTextClassifier":
        """Load a classifier saved with :meth:`save`."""

        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(
            label_counts={key: int(value) for key, value in payload["label_counts"].items()},
            token_counts_by_label={
                label: {token: int(count) for token, count in counts.items()}
                for label, counts in payload["token_counts_by_label"].items()
            },
            token_totals_by_label={
                key: int(value) for key, value in payload["token_totals_by_label"].items()
            },
            vocabulary=set(payload["vocabulary"]),
            alpha=float(payload.get("alpha", 1.0)),
        )
