"""Tests for deepretro.report_scores."""

from __future__ import annotations

import json
from pathlib import Path

from deepretro import report_scores
from deepretro.batch import slugify_molecule

MOL_A = "CCO"
MOL_B = "CC(=O)O"


def _write_pathway(
    base: Path, mode: str, model: str, smi: str, *, solved: bool, sa: float, sc: float
) -> None:
    d = base / mode / model / slugify_molecule(smi)
    d.mkdir(parents=True, exist_ok=True)
    (d / "pathway_1.json").write_text(
        json.dumps(
            {
                "solved": solved,
                "target": smi,
                "scores": {
                    "route_summary": {"mean_product_sa": sa, "mean_product_sc": sc}
                },
            }
        )
    )


def _make_base(tmp_path: Path) -> Path:
    base = tmp_path / "batch"
    _write_pathway(base, "heuristic", "opus", MOL_A, solved=True, sa=2.0, sc=1.5)
    _write_pathway(base, "ml", "opus", MOL_A, solved=False, sa=3.0, sc=2.5)
    # MOL_B: heuristic solved, ml OOM-skipped
    _write_pathway(base, "heuristic", "opus", MOL_B, solved=True, sa=4.0, sc=3.0)
    oom = base / "ml" / "opus" / slugify_molecule(MOL_B)
    oom.mkdir(parents=True, exist_ok=True)
    (oom / "oom_skipped.json").write_text(json.dumps({"smiles": MOL_B}))
    # a hallucination_model training dir that must NOT be treated as a config
    (base / "hallucination_model" / "xgboost_model").mkdir(parents=True, exist_ok=True)
    (base / "hallucination_model" / "xgboost_model" / "model.joblib").write_text("x")
    return base


def test_discover_configs_ignores_training_dir(tmp_path: Path) -> None:
    base = _make_base(tmp_path)
    modes, models = report_scores.discover_configs(base)
    assert modes == ["heuristic", "ml"]  # heuristic ordered first
    assert models == ["opus"]  # xgboost_model excluded


def test_build_report_counts_and_per_mode_scores(tmp_path: Path) -> None:
    base = _make_base(tmp_path)
    t1, t2, labels = report_scores.build_report(
        base, [MOL_A, MOL_B], ["opus"], ["heuristic", "ml"]
    )
    row = t1[0]
    assert row["counts"]["heuristic"] == "2/2"  # both solved
    assert (
        row["counts"]["ml"] == "2/0"
    )  # 1 pathway (unsolved) + 1 oom = 2 done, 0 solved
    # per-mode SA averages are separate
    assert row["sa"]["heuristic"].startswith("3.00")  # mean(2.0, 4.0)
    assert row["sa"]["ml"].startswith("3.00")  # only MOL_A scored (3.0)
    assert labels[MOL_A] == "m01"


def test_generate_writes_md_and_tex(tmp_path: Path) -> None:
    base = _make_base(tmp_path)
    mol_file = tmp_path / "mols.txt"
    mol_file.write_text(f"{MOL_A}\n{MOL_B}\n")
    report_scores.main(
        [
            "--base",
            str(base),
            "--molecules",
            str(mol_file),
            "--models",
            "opus",
            "--modes",
            "heuristic,ml",
        ]
    )
    md = (base / "scores_report.md").read_text()
    tex = (base / "scores_tables.tex").read_text()
    assert "SAScore" in md and "SCScore" in md
    assert "[1, 10]" in md and "[1, 5]" in md  # ranges present
    assert r"\begin{tabular}" in tex
