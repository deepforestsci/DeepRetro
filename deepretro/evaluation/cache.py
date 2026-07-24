"""Small local SQLite cache for reproducible official score reuse."""

from __future__ import annotations

import sqlite3
from pathlib import Path

from deepretro.evaluation.types import ChemCensorProvenance, ChemCensorScore


class ScoreCache:
    """Persist successful official scores under their full scorer provenance.

    Processing failures are intentionally not cached: an unavailable mapper or
    temporary SQLite problem must not turn into a permanent zero-coverage
    artifact in future evaluation runs.
    """

    def __init__(self, path: str | Path) -> None:
        """Create the local cache database and its provenance-aware table."""
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS chemcensor_scores (
                    canonical_reaction_smiles TEXT NOT NULL,
                    scorer_name TEXT NOT NULL,
                    package_version TEXT NOT NULL,
                    database_version TEXT NOT NULL,
                    database_sha256 TEXT NOT NULL,
                    max_center_type INTEGER NOT NULL,
                    find_exact_match INTEGER NOT NULL,
                    with_functional_groups REAL NOT NULL,
                    without_functional_groups REAL NOT NULL,
                    status TEXT NOT NULL,
                    PRIMARY KEY (
                        canonical_reaction_smiles,
                        scorer_name,
                        package_version,
                        database_version,
                        database_sha256,
                        max_center_type,
                        find_exact_match
                    )
                )
                """
            )

    def get(
        self,
        canonical_reaction_smiles: str,
        provenance: ChemCensorProvenance,
        *,
        prediction_id: str,
    ) -> ChemCensorScore | None:
        """Return a matching cached score rebound to the current prediction ID."""
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT with_functional_groups, without_functional_groups, status
                FROM chemcensor_scores
                WHERE canonical_reaction_smiles = ?
                  AND scorer_name = ?
                  AND package_version = ?
                  AND database_version = ?
                  AND database_sha256 = ?
                  AND max_center_type = ?
                  AND find_exact_match = ?
                """,
                _cache_key(canonical_reaction_smiles, provenance),
            ).fetchone()
        if row is None:
            return None
        return ChemCensorScore(
            prediction_id=prediction_id,
            canonical_reaction_smiles=canonical_reaction_smiles,
            with_functional_groups=float(row[0]),
            without_functional_groups=float(row[1]),
            status=str(row[2]),
            provenance=provenance,
        )

    def put(self, score: ChemCensorScore) -> None:
        """Store only stable, successful/no-precedent official results."""
        if score.status not in {"scored", "no_precedent"}:
            return
        if score.provenance is None or score.canonical_reaction_smiles is None:
            raise ValueError("cacheable ChemCensor scores require provenance and reaction")
        assert score.with_functional_groups is not None
        assert score.without_functional_groups is not None
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO chemcensor_scores (
                    canonical_reaction_smiles,
                    scorer_name,
                    package_version,
                    database_version,
                    database_sha256,
                    max_center_type,
                    find_exact_match,
                    with_functional_groups,
                    without_functional_groups,
                    status
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(
                    canonical_reaction_smiles,
                    scorer_name,
                    package_version,
                    database_version,
                    database_sha256,
                    max_center_type,
                    find_exact_match
                ) DO UPDATE SET
                    with_functional_groups = excluded.with_functional_groups,
                    without_functional_groups = excluded.without_functional_groups,
                    status = excluded.status
                """,
                (
                    *_cache_key(score.canonical_reaction_smiles, score.provenance),
                    score.with_functional_groups,
                    score.without_functional_groups,
                    score.status,
                ),
            )

    def _connect(self) -> sqlite3.Connection:
        """Open a short-lived transaction-safe SQLite connection."""
        connection = sqlite3.connect(self.path, timeout=30)
        connection.execute("PRAGMA journal_mode=WAL")
        return connection


def _cache_key(
    canonical_reaction_smiles: str,
    provenance: ChemCensorProvenance,
) -> tuple[str, str, str, str, str, int, int]:
    """Return the exact tuple used by both cache reads and writes."""
    return (
        canonical_reaction_smiles,
        provenance.scorer_name,
        provenance.package_version,
        provenance.database_version,
        provenance.database_sha256,
        provenance.max_center_type,
        int(provenance.find_exact_match),
    )
