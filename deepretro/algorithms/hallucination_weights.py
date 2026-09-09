"""Validated, immutable configuration for structural heuristic scoring."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path


@dataclass(frozen=True)
class HallucinationWeights:
    """Penalty weights, caps, and severity cutoffs for the heuristic checker.

    Parameters
    ----------
    w_atom : int
        Points deducted per mismatched atom, per element.
    cap_atom : int
        Maximum atom-mismatch deduction.
    w_ring : int
        Points deducted per ring-size change.
    cap_ring : int
        Maximum ring-change deduction.
    w_pos : int
        Points deducted per substituent position change.
    cap_pos : int
        Maximum substituent-position deduction.
    w_arom : int
        Flat deduction when aromaticity changes significantly.
    arom_trigger : int
        The aromaticity penalty fires when the absolute change in aromatic atom
        count is strictly greater than this. A threshold on an atom count, not a
        penalty, so it is never rescaled.
    w_bond : int
        Points deducted per unexplained extra bond.
    cap_bond : int
        Maximum extra-bond deduction.
    cut_low : int
        Scores at or above this are severity ``"low"``.
    cut_medium : int
        Scores at or above this (and below ``cut_low``) are ``"medium"``. The
        display cutoff is independent of candidate acceptance.
    cut_high : int
        Scores at or above this (and below ``cut_medium``) are ``"high"``.
        Anything lower is ``"critical"``.

    reject_below : int
        Flag scores below this threshold. Ranking retains flagged candidates
        as fallback; hard filtering excludes them. Range 0 through 101.

    Raises
    ------
    TypeError
        If any field is not an ``int``.
    ValueError
        If any field is negative, or if the cutoffs are not strictly ordered
        ``cut_high < cut_medium < cut_low``.

    Examples
    --------
    >>> from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS
    >>> DEFAULT_WEIGHTS.w_atom
    5
    >>> DEFAULT_WEIGHTS.replace(w_atom=9).w_atom
    9
    """

    w_atom: int = 5
    cap_atom: int = 100
    w_ring: int = 25
    cap_ring: int = 50
    w_pos: int = 60
    cap_pos: int = 100
    w_arom: int = 40
    arom_trigger: int = 2
    w_bond: int = 5
    cap_bond: int = 30
    cut_low: int = 80
    cut_medium: int = 40
    cut_high: int = 20
    #: Acceptance threshold, independent of display severity.
    reject_below: int = 40

    def __post_init__(self) -> None:
        """Validate types, non-negativity, and cutoff ordering."""
        for field in fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, int) or isinstance(value, bool):
                raise TypeError(f"{field.name} must be an int, got {value!r}")
            if value < 0:
                raise ValueError(f"{field.name} must be non-negative, got {value}")
        if self.reject_below > 101:
            raise ValueError(
                f"reject_below must be <= 101 (scores are 0-100), got "
                f"{self.reject_below}"
            )
        if not (self.cut_high < self.cut_medium < self.cut_low):
            raise ValueError(
                "cutoffs must satisfy cut_high < cut_medium < cut_low, got "
                f"{self.cut_high} < {self.cut_medium} < {self.cut_low}"
            )

    def replace(self, **kwargs: int) -> HallucinationWeights:
        """Return a new validated instance with fields overridden.

        Parameters
        ----------
        **kwargs : int
            Field names and their new values.

        Returns
        -------
        HallucinationWeights
            A new instance.

        Examples
        --------
        >>> from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS
        >>> DEFAULT_WEIGHTS.replace(w_ring=30).w_ring
        30
        """
        return replace(self, **kwargs)

    def to_json(self, path: str | Path) -> None:
        """Write all weight fields to ``path`` as JSON.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination; parent directories are created.

        Returns
        -------
        None

        Examples
        --------
        >>> import tempfile, pathlib
        >>> from deepretro.algorithms.hallucination_weights import DEFAULT_WEIGHTS
        >>> p = pathlib.Path(tempfile.mkdtemp()) / "w.json"
        >>> DEFAULT_WEIGHTS.to_json(p)
        >>> p.exists()
        True
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(asdict(self), indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    @classmethod
    def from_json(cls, path: str | Path) -> HallucinationWeights:
        """Load an instance written by :meth:`to_json`.

        The file must hold all weight fields. Legacy files without reject_below
        inherit their cut_medium threshold. Unknown fields are
        rejected so tuning metadata cannot leak in, and missing fields are
        rejected so a truncated file cannot silently fall back to defaults.

        Parameters
        ----------
        path : str or pathlib.Path
            Source file.

        Returns
        -------
        HallucinationWeights
            The validated instance.

        Raises
        ------
        ValueError
            If the file carries fields that are not weights, or omits any.

        Examples
        --------
        >>> import tempfile, pathlib
        >>> from deepretro.algorithms.hallucination_weights import (
        ...     DEFAULT_WEIGHTS, HallucinationWeights)
        >>> p = pathlib.Path(tempfile.mkdtemp()) / "w.json"
        >>> DEFAULT_WEIGHTS.to_json(p)
        >>> HallucinationWeights.from_json(p) == DEFAULT_WEIGHTS
        True
        """
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("weight JSON must contain an object")
        known = {field.name for field in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"unknown weight fields in {path}: {sorted(unknown)}")

        # reject_below post-dates the first saved weight files. Those gated on
        # cut_medium, so defaulting to it reproduces the file's original meaning
        # exactly rather than silently applying today's default of 40.
        if "reject_below" not in data and "cut_medium" in data:
            data = {**data, "reject_below": data["cut_medium"]}

        missing = known - set(data)
        if missing:
            raise ValueError(f"missing weight fields in {path}: {sorted(missing)}")
        return cls(**data)


DEFAULT_WEIGHTS = HallucinationWeights()
