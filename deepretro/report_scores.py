"""Summarise SA / SCScore across a batch run into Markdown and LaTeX tables.

Reads the pathway JSONs written by :mod:`deepretro.batch` (layout
``<base>/<mode>/<model>/<molecule_slug>/pathway_*.json``) and emits two tables:

1. Per model: heuristic and ml done/solved counts, plus average (max, min)
   SAScore and SCScore.
2. Per molecule x model: the ml and heuristic SA/SC for each route.

Both scores are "lower is better" (SA: more synthetically accessible; SC: less
complex), rendered with a downward arrow. The tool only *reads* scores already
present in the pathway JSONs (produced by :func:`deepretro.score.score_pathway`);
it does not recompute them.

Example
-------
    python scripts/generate_score_report.py \\
        --base ~/Downloads/deepretro_batch_2026-07-22_150635 \\
        --molecules molecules.txt
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import mean
from typing import Any

from deepretro.batch import read_molecules, slugify_molecule

_TEX_ESCAPE = {
    "\\": r"\textbackslash{}",
    "_": r"\_",
    "#": r"\#",
    "$": r"\$",
    "%": r"\%",
    "&": r"\&",
    "{": r"\{",
    "}": r"\}",
    "^": r"\^{}",
    "~": r"\~{}",
}


def _tex_escape(text: str) -> str:
    """Escape LaTeX-special characters in a string."""
    return "".join(_TEX_ESCAPE.get(ch, ch) for ch in text)


def _fmt(value: Any) -> str:
    """Format a number to 2 dp, or ``-`` when unavailable."""
    return f"{value:.2f}" if isinstance(value, (int, float)) else "-"


def discover_configs(base: Path) -> tuple[list[str], list[str]]:
    """Discover the ``(modes, models)`` present under *base*.

    Parameters
    ----------
    base : Path
        Batch output root (``<base>/<mode>/<model>/...``).

    Returns
    -------
    tuple[list[str], list[str]]
        Sorted modes (``heuristic`` first if present) and sorted models.

    Examples
    --------
    >>> discover_configs(Path("/nonexistent"))
    ([], [])
    """
    if not base.is_dir():
        return [], []

    def _has_outputs(model_dir: Path) -> bool:
        # A real config dir contains molecule subdirs with pathway/oom/error
        # files (this excludes e.g. the hallucination_model/ training dir).
        for mol in model_dir.iterdir():
            if not mol.is_dir():
                continue
            if (
                any(mol.glob("pathway_*.json"))
                or (mol / "oom_skipped.json").exists()
                or (mol / "error.json").exists()
            ):
                return True
        return False

    modes: set[str] = set()
    models: set[str] = set()
    for mode_dir in (p for p in base.iterdir() if p.is_dir()):
        for model_dir in (p for p in mode_dir.iterdir() if p.is_dir()):
            if _has_outputs(model_dir):
                modes.add(mode_dir.name)
                models.add(model_dir.name)
    ordered_modes = sorted(modes, key=lambda m: (m != "heuristic", m))
    return ordered_modes, sorted(models)


def load_cell(base: Path, mode: str, model: str, smiles: str) -> dict[str, Any]:
    """Load one (mode, model, molecule) outcome from disk.

    Returns a dict with ``state`` in ``{"done", "oom", "err", "todo"}`` and,
    when done, ``solved`` plus ``sa``/``sc`` route means.
    """
    mol_dir = base / mode / model / slugify_molecule(smiles)
    pathway = mol_dir / "pathway_1.json"
    if pathway.exists():
        payload = json.loads(pathway.read_text(encoding="utf-8"))
        summary = payload.get("scores", {}).get("route_summary", {})
        return {
            "state": "done",
            "solved": bool(payload.get("solved")),
            "sa": summary.get("mean_product_sa"),
            "sc": summary.get("mean_product_sc"),
        }
    if (mol_dir / "oom_skipped.json").exists():
        return {"state": "oom"}
    if (mol_dir / "error.json").exists():
        return {"state": "err"}
    return {"state": "todo"}


def _cell_value(cell: dict[str, Any], key: str) -> str:
    if cell["state"] == "done":
        return _fmt(cell.get(key))
    return {"oom": "OOM", "err": "ERR", "todo": "-"}[cell["state"]]


def build_report(
    base: Path,
    molecules: list[str],
    models: list[str],
    modes: list[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, str]]:
    """Build the per-model (table 1) and per-molecule (table 2) rows.

    Returns
    -------
    tuple
        ``(table1_rows, table2_rows, molecule_labels)``.
    """
    labels = {smi: f"m{i + 1:02d}" for i, smi in enumerate(molecules)}
    data = {
        (model, mode, smi): load_cell(base, mode, model, smi)
        for model in models
        for mode in modes
        for smi in molecules
    }

    def _agg(values: list[float]) -> str:
        return (
            f"{_fmt(mean(values))} ({_fmt(max(values))}, {_fmt(min(values))})"
            if values
            else "-"
        )

    table1: list[dict[str, Any]] = []
    for model in models:
        counts: dict[str, str] = {}
        sa_by_mode: dict[str, str] = {}
        sc_by_mode: dict[str, str] = {}
        for mode in modes:
            cells = [data[(model, mode, s)] for s in molecules]
            done = sum(1 for c in cells if c["state"] in ("done", "oom", "err"))
            solved = sum(1 for c in cells if c.get("solved"))
            counts[mode] = f"{done}/{solved}"
            sa_by_mode[mode] = _agg([c["sa"] for c in cells if c.get("sa") is not None])
            sc_by_mode[mode] = _agg([c["sc"] for c in cells if c.get("sc") is not None])
        table1.append(
            {
                "model": model,
                "counts": counts,
                "sa": sa_by_mode,
                "sc": sc_by_mode,
            }
        )

    table2: list[dict[str, Any]] = []
    for smi in molecules:
        for model in models:
            row = {"mol": labels[smi], "model": model}
            for mode in modes:
                cell = data[(model, mode, smi)]
                row[f"{mode}_sa"] = _cell_value(cell, "sa")
                row[f"{mode}_sc"] = _cell_value(cell, "sc")
            table2.append(row)

    return table1, table2, labels


def render_markdown(
    base: Path,
    table1: list[dict[str, Any]],
    table2: list[dict[str, Any]],
    labels: dict[str, str],
    modes: list[str],
) -> str:
    """Render the report as Markdown."""
    lines = [
        "# DeepRetro batch SA / SCScore report",
        f"\n_Source: `{base}`_\n",
        "**Score ranges** (both **lower is better**):",
        "- **SAScore** ∈ [1, 10] — 1 = easy to synthesize, 10 = very hard (↓ better).",
        "- **SCScore** ∈ [1, 5] — 1 = simple, 5 = complex (↓ better).\n",
        "## Table 1 — per model",
        "| Model | "
        + " | ".join(f"{m} done/solved" for m in modes)
        + " | "
        + " | ".join(f"{m} SA avg (max, min) ↓" for m in modes)
        + " | "
        + " | ".join(f"{m} SC avg (max, min) ↓" for m in modes)
        + " |",
        "|---" * (1 + 3 * len(modes)) + "|",
    ]
    for r in table1:
        counts = " | ".join(r["counts"].get(m, "-") for m in modes)
        sa = " | ".join(r["sa"].get(m, "-") for m in modes)
        sc = " | ".join(r["sc"].get(m, "-") for m in modes)
        lines.append(f"| {r['model']} | {counts} | {sa} | {sc} |")

    lines.append("\n## Table 2 — per molecule × model")
    head = (
        "| Molecule | Model | "
        + " | ".join(f"{m}-SA ↓ | {m}-SC ↓" for m in modes)
        + " |"
    )
    lines.append(head)
    lines.append("|---" * (2 + 2 * len(modes)) + "|")
    for r in table2:
        vals = " | ".join(f"{r[f'{m}_sa']} | {r[f'{m}_sc']}" for m in modes)
        lines.append(f"| {r['mol']} | {r['model']} | {vals} |")

    lines.append("\n## Molecule legend")
    for smi, label in labels.items():
        lines.append(f"- **{label}** — `{smi}`")
    return "\n".join(lines) + "\n"


def render_latex(
    table1: list[dict[str, Any]],
    table2: list[dict[str, Any]],
    labels: dict[str, str],
    modes: list[str],
) -> str:
    """Render the two tables as standalone LaTeX ``table`` environments."""
    t1_cols = "l" + "c" * (3 * len(modes))
    t1_head = (
        "Model & "
        + " & ".join(f"{_tex_escape(m)} (d/s)" for m in modes)
        + " & "
        + " & ".join(rf"{_tex_escape(m)} SA (max,min)$\downarrow$" for m in modes)
        + " & "
        + " & ".join(rf"{_tex_escape(m)} SC (max,min)$\downarrow$" for m in modes)
        + r" \\"
    )
    lines = [
        "% DeepRetro SA/SCScore tables. Both lower is better.",
        "% SAScore range [1,10] (1 easy, 10 hard); SCScore range [1,5] "
        "(1 simple, 5 complex).",
        r"\begin{table}[htbp]\centering\small",
        r"\caption{Per-model summary; SA$\downarrow$ and SC$\downarrow$ split by "
        r"hallucination mode.}",
        rf"\begin{{tabular}}{{{t1_cols}}}",
        r"\hline",
        t1_head,
        r"\hline",
    ]
    for r in table1:
        counts = " & ".join(r["counts"].get(m, "-") for m in modes)
        sa = " & ".join(_tex_escape(r["sa"].get(m, "-")) for m in modes)
        sc = " & ".join(_tex_escape(r["sc"].get(m, "-")) for m in modes)
        lines.append(f"{_tex_escape(r['model'])} & {counts} & {sa} & {sc} \\\\")
    lines += [r"\hline", r"\end{tabular}\end{table}", ""]

    t2_cols = "ll" + "c" * (2 * len(modes))
    t2_head = (
        "Molecule & Model & "
        + " & ".join(f"{_tex_escape(m)}-SA & {_tex_escape(m)}-SC" for m in modes)
        + r" \\"
    )
    lines += [
        r"\begin{table}[htbp]\centering\small",
        r"\caption{Per-molecule $\times$ model SA/SC ($\downarrow$ better).}",
        rf"\begin{{tabular}}{{{t2_cols}}}",
        r"\hline",
        t2_head,
        r"\hline",
    ]
    for r in table2:
        vals = " & ".join(f"{r[f'{m}_sa']} & {r[f'{m}_sc']}" for m in modes)
        lines.append(f"{r['mol']} & {_tex_escape(r['model'])} & {vals} \\\\")
    lines += [r"\hline", r"\end{tabular}\end{table}", "", "% Molecule legend:"]
    for smi, label in labels.items():
        lines.append(f"% {label} = {smi}")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> None:
    """Generate the Markdown + LaTeX score report from the command line."""
    parser = argparse.ArgumentParser(description="DeepRetro SA/SCScore report")
    parser.add_argument("--base", required=True, help="Batch output root directory")
    parser.add_argument(
        "--molecules", required=True, help="Molecules file (one SMILES/line)"
    )
    parser.add_argument(
        "--models", default=None, help="Comma-separated models (default: auto-discover)"
    )
    parser.add_argument(
        "--modes", default=None, help="Comma-separated modes (default: auto-discover)"
    )
    parser.add_argument(
        "--out-md",
        default=None,
        help="Markdown output path (default: <base>/scores_report.md)",
    )
    parser.add_argument(
        "--out-tex",
        default=None,
        help="LaTeX output path (default: <base>/scores_tables.tex)",
    )
    args = parser.parse_args(argv)

    base = Path(args.base).expanduser()
    disc_modes, disc_models = discover_configs(base)
    modes = args.modes.split(",") if args.modes else (disc_modes or ["heuristic", "ml"])
    models = args.models.split(",") if args.models else disc_models
    molecules = read_molecules(args.molecules)

    table1, table2, labels = build_report(base, molecules, models, modes)
    out_md = Path(args.out_md) if args.out_md else base / "scores_report.md"
    out_tex = Path(args.out_tex) if args.out_tex else base / "scores_tables.tex"
    out_md.write_text(
        render_markdown(base, table1, table2, labels, modes), encoding="utf-8"
    )
    out_tex.write_text(render_latex(table1, table2, labels, modes), encoding="utf-8")
    print(f"Wrote {out_md}")
    print(f"Wrote {out_tex}")


if __name__ == "__main__":  # pragma: no cover
    main()
