"""Tests for the batch retrosynthesis runner."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from deepretro import batch


class FakeResponse:
    """Minimal stand-in for a requests.Response."""

    def __init__(self, content: bytes, status: int = 200) -> None:
        self.content = content
        self.status = status

    def raise_for_status(self) -> None:
        if self.status >= 400:
            raise RuntimeError(f"HTTP {self.status}")


class FakeSolver:
    """AutoSolver-shaped double returning two routes for any molecule."""

    def solve_multiple(self, smiles: str, k: int = 3) -> list[tuple[dict, bool]]:
        return [
            ({"type": "mol", "smiles": smiles, "children": []}, True),
            ({"type": "mol", "smiles": smiles, "children": []}, False),
        ][:k]

    def parse(self, route: dict, *, solved: bool) -> dict[str, Any]:
        return {"steps": [], "dependencies": {}, "solved": solved}

    def add_metadata(self, parsed: dict[str, Any]) -> dict[str, Any]:
        parsed["enriched"] = True
        return parsed


def test_read_molecules_skips_blanks_and_comments(tmp_path: Path) -> None:
    path = tmp_path / "mols.txt"
    path.write_text("CCO\n\n# a comment\n  CC(=O)O  \n")
    assert batch.read_molecules(str(path)) == ["CCO", "CC(=O)O"]


def test_slugify_is_filesystem_safe_and_deterministic() -> None:
    slug = batch.slugify_molecule("CC(=O)Oc1ccccc1C(=O)O")
    assert "/" not in slug and "(" not in slug and ")" not in slug
    assert slug == batch.slugify_molecule("CC(=O)Oc1ccccc1C(=O)O")


def test_slugify_distinguishes_different_molecules() -> None:
    assert batch.slugify_molecule("CCO") != batch.slugify_molecule("CCN")


def test_download_sheet_csv_writes_content(tmp_path: Path) -> None:
    dest = tmp_path / "data.csv"

    def fake_get(url: str, timeout: float | None = None) -> FakeResponse:
        return FakeResponse(b"product,reactants,label\nCCO,CC,0\n")

    batch.download_sheet_csv("http://sheet", str(dest), http_get=fake_get)
    assert dest.read_bytes().startswith(b"product,reactants,label")


def test_download_sheet_csv_raises_on_http_error(tmp_path: Path) -> None:
    def fake_get(url: str, timeout: float | None = None) -> FakeResponse:
        return FakeResponse(b"", status=404)

    with pytest.raises(RuntimeError):
        batch.download_sheet_csv(
            "http://sheet", str(tmp_path / "x.csv"), http_get=fake_get
        )


def test_solve_molecule_attaches_target_to_each_pathway() -> None:
    pathways = batch.solve_molecule(FakeSolver(), "CCO", top_k=3)
    assert len(pathways) == 2
    assert all(p["target"] == "CCO" for p in pathways)
    assert all(p["enriched"] is True for p in pathways)


def test_run_batch_writes_pathway_files(tmp_path: Path) -> None:
    def solve(smiles: str) -> list[dict[str, Any]]:
        return [{"target": smiles, "solved": True}, {"target": smiles, "solved": False}]

    batch.run_batch(
        ["CCO"], str(tmp_path), timestamp="2026-07-01_00-00-00", solve=solve
    )
    mol_dir = tmp_path / "2026-07-01_00-00-00" / batch.slugify_molecule("CCO")
    assert (mol_dir / "pathway_1.json").exists()
    assert (mol_dir / "pathway_2.json").exists()
    payload = json.loads((mol_dir / "pathway_1.json").read_text())
    assert payload["target"] == "CCO"


def test_run_batch_writes_error_json_on_failure(tmp_path: Path) -> None:
    def solve(smiles: str) -> list[dict[str, Any]]:
        raise RuntimeError("solver exploded")

    batch.run_batch(
        ["CCO"], str(tmp_path), timestamp="2026-07-01_00-00-00", solve=solve
    )
    error_file = (
        tmp_path / "2026-07-01_00-00-00" / batch.slugify_molecule("CCO") / "error.json"
    )
    assert error_file.exists()
    assert "solver exploded" in error_file.read_text()


def test_train_hallucination_checker_skips_without_label_column(tmp_path: Path) -> None:
    csv_path = tmp_path / "unlabeled.csv"
    csv_path.write_text("product,reactants\nCCO,CC\n")
    result = batch.train_hallucination_checker(str(csv_path), str(tmp_path / "model"))
    assert result is None
