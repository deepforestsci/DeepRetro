"""Tests for the batch retrosynthesis runner."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from deepretro import batch
from deepretro.utils.llm_trace import LOG_FILENAME, current_trace, record_llm_call


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


def test_arg_parser_exposes_agent_iteration_budget_flags() -> None:
    args = batch._build_arg_parser().parse_args(["--molecules", "m.txt"])
    assert args.agent_min_iterations == 5
    assert args.agent_max_iterations == 15
    assert args.agent_iteration_decay == 0.75


def test_main_passes_agent_iteration_budget_to_solver(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import deepretro.algorithms.autosolve as autosolve_module

    captured: dict[str, Any] = {}

    class RecordingSolver:
        def __init__(self, **kwargs: Any) -> None:
            captured.update(kwargs)

    monkeypatch.setattr(autosolve_module, "AutoSolver", RecordingSolver)
    monkeypatch.setattr(batch, "run_batch", lambda *a, **kw: None)
    molecules = tmp_path / "m.txt"
    molecules.write_text("CCO\n")

    batch.main(
        [
            "--molecules",
            str(molecules),
            "--out",
            str(tmp_path / "out"),
            "--hallucination-mode",
            "none",
            "--agent-min-iterations",
            "4",
            "--agent-max-iterations",
            "12",
            "--agent-iteration-decay",
            "0.5",
        ]
    )

    assert captured["agent_min_iterations"] == 4
    assert captured["agent_max_iterations"] == 12
    assert captured["agent_iteration_decay"] == 0.5


def test_run_batch_writes_llm_call_log_in_the_molecule_dir(tmp_path: Path) -> None:
    """LLM calls made while solving land in the molecule's ``llm_calls.jsonl``."""

    def solve(smiles: str) -> list[dict[str, Any]]:
        record_llm_call(
            stage="retrosynthesis",
            model="openai/gpt-4o-mini",
            messages=[{"role": "user", "content": smiles}],
            response="ok",
            latency_ms=2.0,
        )
        return [{"target": smiles, "solved": True}]

    batch.run_batch(
        ["CCO"], str(tmp_path), timestamp="2026-07-01_00-00-00", solve=solve
    )
    log_path = (
        tmp_path / "2026-07-01_00-00-00" / batch.slugify_molecule("CCO") / LOG_FILENAME
    )
    lines = log_path.read_text().splitlines()
    assert len(lines) == 1
    record = json.loads(lines[0])
    assert record["target"] == "CCO"
    assert record["stage"] == "retrosynthesis"


def test_run_batch_trace_is_active_when_the_solver_raises(tmp_path: Path) -> None:
    """A failing molecule still writes its call log and its error.json."""

    def solve(smiles: str) -> list[dict[str, Any]]:
        record_llm_call(
            stage="retrosynthesis",
            model="m",
            messages=[],
            response=None,
            error="provider down",
            latency_ms=1.0,
        )
        raise RuntimeError("solver exploded")

    batch.run_batch(
        ["CCO"], str(tmp_path), timestamp="2026-07-01_00-00-00", solve=solve
    )
    mol_dir = tmp_path / "2026-07-01_00-00-00" / batch.slugify_molecule("CCO")
    assert "solver exploded" in (mol_dir / "error.json").read_text()
    record = json.loads((mol_dir / LOG_FILENAME).read_text().splitlines()[0])
    assert record["error"] == "provider down"


def test_run_batch_gives_each_molecule_its_own_session(tmp_path: Path) -> None:
    """Two targets never share a Langfuse session id or a log file."""
    sessions: list[str] = []

    def solve(smiles: str) -> list[dict[str, Any]]:
        trace = current_trace()
        assert trace is not None
        sessions.append(trace.session_id)
        return [{"target": smiles}]

    batch.run_batch(
        ["CCO", "CCC"], str(tmp_path), timestamp="2026-07-01_00-00-00", solve=solve
    )
    assert len(set(sessions)) == 2
