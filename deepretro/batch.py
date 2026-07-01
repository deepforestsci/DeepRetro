"""Batch retrosynthesis runner.

End-to-end driver: download a CSV from a public Google Sheet, (optionally) train
the hallucination checker, read target molecules from a text file, run each
through :class:`deepretro.algorithms.autosolve.AutoSolver`, and dump the routes
as ``<out>/<timestamp>/<molecule>/pathway_<i>.json``.

The training step is a **template** (see :func:`train_hallucination_checker`): it
runs only when the CSV carries the expected labelled columns, and otherwise logs
and falls back to the heuristic hallucination mode so the batch still runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import structlog

logger = structlog.get_logger(__name__)

# Solver hook: molecule SMILES -> list of parsed pathway dicts.
SolveMolecule = Callable[[str], list[dict[str, Any]]]
# HTTP hook: (url, timeout) -> response with ``content`` and ``raise_for_status``.
HttpGet = Callable[..., Any]

REQUIRED_TRAINING_COLUMNS = {"product", "reactants", "label"}


def read_molecules(path: str) -> list[str]:
    """Read target SMILES from a text file (one per line).

    Blank lines and lines starting with ``#`` are skipped; surrounding
    whitespace is stripped.

    Parameters
    ----------
    path : str
        Path to the molecules file.

    Returns
    -------
    list[str]
        Target SMILES strings.

    Examples
    --------
    >>> import tempfile, os
    >>> path = os.path.join(tempfile.mkdtemp(), "m.txt")
    >>> _ = open(path, "w").write("CCO\\n# note\\n\\nCCN\\n")
    >>> read_molecules(path)
    ['CCO', 'CCN']
    """
    molecules: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        molecules.append(stripped)
    return molecules


def slugify_molecule(smiles: str) -> str:
    """Build a filesystem-safe, unique directory name for a molecule.

    Combines a readable ASCII-safe prefix with a short hash of the original
    SMILES so distinct molecules never collide and the same SMILES is stable.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.

    Returns
    -------
    str
        A safe directory name.

    Examples
    --------
    >>> slugify_molecule("CCO") == slugify_molecule("CCO")
    True
    >>> "/" in slugify_molecule("CC(=O)O")
    False
    """
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", smiles)[:60]
    digest = hashlib.md5(smiles.encode("utf-8")).hexdigest()[:8]
    return f"{safe}_{digest}"


def download_sheet_csv(
    url: str,
    dest: str,
    *,
    http_get: HttpGet | None = None,
) -> str:
    """Download a published Google Sheet CSV export to *dest*.

    Parameters
    ----------
    url : str
        The sheet's ``/export?format=csv`` URL (published / link-shared; no
        auth is used).
    dest : str
        Local path to write the CSV to.
    http_get : callable, optional
        Injectable HTTP getter (defaults to ``requests.get``). Must return a
        response exposing ``content`` and ``raise_for_status``.

    Returns
    -------
    str
        The destination path.

    Raises
    ------
    Exception
        Whatever ``raise_for_status`` raises on a non-2xx response.
    """
    if http_get is None:
        import requests

        http_get = requests.get
    response = http_get(url, timeout=60)
    response.raise_for_status()
    Path(dest).write_bytes(response.content)
    logger.info("Downloaded sheet CSV", url=url, dest=dest)
    return dest


def train_hallucination_checker(csv_path: str, save_dir: str) -> str | None:
    """Train the ML hallucination checker from a labelled CSV (TEMPLATE).

    Runs only when the CSV has ``product``/``reactants``/``label`` columns;
    otherwise logs and returns ``None`` so the batch falls back to the
    heuristic checker. The training body is a template for a downstream
    contributor to complete/tune for their dataset; any failure degrades to
    ``None`` rather than aborting the batch.

    Parameters
    ----------
    csv_path : str
        Path to the training CSV.
    save_dir : str
        Directory to save the trained model into.

    Returns
    -------
    str or None
        ``save_dir`` when a model was trained and saved, else ``None``.
    """
    import pandas as pd

    columns = set(pd.read_csv(csv_path, nrows=0).columns)
    if not REQUIRED_TRAINING_COLUMNS.issubset(columns):
        logger.warning(
            "Training CSV missing required columns; skipping training (template)",
            required=sorted(REQUIRED_TRAINING_COLUMNS),
            found=sorted(columns),
        )
        return None

    try:
        # --- TEMPLATE: adapt featurization/splitting/hyperparameters as needed ---
        from deepretro.data.loader import ReactionDataLoader, stratified_split
        from deepretro.models import HallucinationClassifier

        loader = ReactionDataLoader()
        dataset = loader.create_dataset(csv_path)
        train_ds, _valid_ds, test_ds = stratified_split(dataset)

        Path(save_dir).mkdir(parents=True, exist_ok=True)
        classifier = HallucinationClassifier(model_dir=save_dir)
        classifier.fit(train_ds)
        classifier.evaluate(test_ds)
        classifier.save(save_dir)
        # --- end template ---
        logger.info("Trained hallucination checker", save_dir=save_dir)
        return save_dir
    except Exception as exc:  # ponytail: template; never abort the batch on training
        logger.error("Hallucination training failed; using heuristic", error=str(exc))
        return None


def solve_molecule(solver: Any, smiles: str, top_k: int) -> list[dict[str, Any]]:
    """Solve one molecule into up to ``top_k`` parsed, enriched pathway dicts.

    Parameters
    ----------
    solver : AutoSolver
        Any object exposing ``solve_multiple``, ``parse``, and ``add_metadata``.
    smiles : str
        Target molecule SMILES.
    top_k : int
        Maximum number of candidate routes.

    Returns
    -------
    list[dict]
        Parsed pathway dicts, each tagged with its ``target`` SMILES.
    """
    pathways: list[dict[str, Any]] = []
    for route, solved in solver.solve_multiple(smiles, k=top_k):
        parsed = solver.parse(route, solved=solved)
        parsed = solver.add_metadata(parsed)
        parsed["target"] = smiles
        pathways.append(parsed)
    return pathways


def run_batch(
    molecules: list[str],
    out_dir: str,
    *,
    timestamp: str,
    solve: SolveMolecule,
) -> dict[str, list[str]]:
    """Run each molecule and write its routes under ``out_dir/timestamp/molecule``.

    A per-molecule failure writes ``error.json`` and never aborts the batch.

    Parameters
    ----------
    molecules : list[str]
        Target SMILES.
    out_dir : str
        Root output directory.
    timestamp : str
        Run timestamp used as the second path level.
    solve : callable
        ``smiles -> list[dict]`` producing the parsed pathways for a molecule.

    Returns
    -------
    dict[str, list[str]]
        Map from each molecule to the file paths written for it.
    """
    run_dir = Path(out_dir) / timestamp
    written: dict[str, list[str]] = {}

    for smiles in molecules:
        mol_dir = run_dir / slugify_molecule(smiles)
        mol_dir.mkdir(parents=True, exist_ok=True)
        try:
            pathways = solve(smiles)
            paths: list[str] = []
            for index, pathway in enumerate(pathways, start=1):
                target = mol_dir / f"pathway_{index}.json"
                target.write_text(json.dumps(pathway, indent=2), encoding="utf-8")
                paths.append(str(target))
            written[smiles] = paths
            logger.info("Wrote pathways", molecule=smiles, count=len(paths))
        except Exception as exc:  # ponytail: one bad molecule must not kill the run
            error_file = mol_dir / "error.json"
            error_file.write_text(
                json.dumps({"smiles": smiles, "error": str(exc)}, indent=2),
                encoding="utf-8",
            )
            written[smiles] = [str(error_file)]
            logger.error("Molecule failed", molecule=smiles, error=str(exc))

    return written


def _build_arg_parser() -> argparse.ArgumentParser:
    """Build the batch runner argument parser."""
    parser = argparse.ArgumentParser(
        description="DeepRetro batch retrosynthesis runner"
    )
    parser.add_argument(
        "--sheet-url",
        required=True,
        help="Google Sheets /export?format=csv URL (public/link-shared). "
        "TODO: supply the sheet URL.",
    )
    parser.add_argument(
        "--molecules",
        required=True,
        help="Path to a text file of target SMILES (one per line). "
        "TODO: supply the molecules file.",
    )
    parser.add_argument("--out", default="batch_output", help="Output directory root")
    parser.add_argument(
        "--csv", default="hallucination_data.csv", help="CSV download path"
    )
    parser.add_argument(
        "--solve-mode",
        choices=["pipeline", "single_step_agent"],
        default="pipeline",
        help="Retrosynthesis mode (orchestrator is not offered)",
    )
    parser.add_argument(
        "--tool-backend", choices=["structured", "sandbox"], default="structured"
    )
    parser.add_argument("--model", default="anthropic/claude-sonnet-4-6")
    parser.add_argument("--az-model", default="USPTO")
    parser.add_argument(
        "--classifier",
        default=None,
        help="Path to a trained hallucination model dir (enables ML mode)",
    )
    parser.add_argument("--top-k", type=int, default=3)
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the batch pipeline from the command line.

    Parameters
    ----------
    argv : list[str], optional
        Argument vector (defaults to ``sys.argv``).
    """
    args = _build_arg_parser().parse_args(argv)
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")

    download_sheet_csv(args.sheet_url, args.csv)

    classifier_dir = args.classifier or train_hallucination_checker(
        args.csv, str(Path(args.out) / "hallucination_model")
    )
    hallucination_mode = "ml" if classifier_dir else "heuristic"

    molecules = read_molecules(args.molecules)
    logger.info(
        "Starting batch",
        molecules=len(molecules),
        timestamp=timestamp,
        solve_mode=args.solve_mode,
        hallucination_mode=hallucination_mode,
    )

    from deepretro.algorithms.autosolve import AutoSolver

    solver = AutoSolver(
        llm=args.model,
        az_model=args.az_model,
        solve_mode=args.solve_mode,
        tool_backend=args.tool_backend,
        hallucination_mode=hallucination_mode,
        hallucination_classifier=classifier_dir,
    )

    run_batch(
        molecules,
        args.out,
        timestamp=timestamp,
        solve=lambda smiles: solve_molecule(solver, smiles, args.top_k),
    )


if __name__ == "__main__":  # pragma: no cover
    main()
