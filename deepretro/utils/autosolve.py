"""Batch and single-molecule retrosynthesis runners.

Wraps the full DeepRetro pipeline (AiZynthFinder -> LLM fallback ->
metadata enrichment) into callable functions for headless execution
without the GUI or server.

* :func:`autosolve` — synchronous, one molecule at a time.
* :func:`autosolve_async` — concurrent execution over a list of SMILES
  using ``asyncio`` + a thread-pool (the underlying pipeline is
  synchronous, so each molecule gets its own worker thread).
"""

from __future__ import annotations

import asyncio
import functools
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import structlog

from src.main import main as _run_retrosynthesis
from src.cache import clear_cache_for_molecule

log = structlog.get_logger(__name__)


def _safe_filename(smiles: str) -> str:
    """Derive a filesystem-safe name from a SMILES string."""
    h = hashlib.md5(smiles.encode()).hexdigest()[:10]
    safe = smiles.replace("/", "_").replace("\\", "_").replace(":", "_")
    if len(safe) > 80:
        safe = safe[:80]
    return f"{safe}_{h}"


def autosolve(
    smiles: str,
    *,
    llm: str = "anthropic/claude-3-7-sonnet-20250219:adv",
    az_model: str = "Pistachio_100+",
    stability_flag: str = "True",
    hallucination_check: str = "True",
    use_protecting_group_feature: bool = False,
    clear_cache: bool = True,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Run retrosynthesis on a single molecule.

    This is the synchronous, single-molecule entry point that replaces
    the old ``prod_dfs_runner.py`` script.  It calls the same pipeline
    — :func:`src.main.main` -> :func:`src.prithvi.run_prithvi` ->
    :func:`src.rec_prithvi.rec_run_prithvi` — but returns the result
    as a dict instead of writing it to a hard-coded path.

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.
    llm : str
        LLM model identifier passed to the pipeline.
    az_model : str
        AiZynthFinder model variant (e.g. ``"USPTO"``,
        ``"Pistachio_100+"``).
    stability_flag : str
        ``"True"`` to enable stability checking.
    hallucination_check : str
        ``"True"`` to enable hallucination checking.
    use_protecting_group_feature : bool
        Enable protecting-group masking in the LLM prompt.
    clear_cache : bool
        Clear any cached AZ / LLM results for this molecule before
        running.  Set ``False`` to reuse previous results.
    output_path : str or Path or None
        If given, the result dict is also written as JSON to this file.

    Returns
    -------
    dict[str, Any]
        Retrosynthesis result dictionary.

    Examples
    --------
    >>> from deepretro.utils.autosolve import autosolve  # doctest: +SKIP
    >>> result = autosolve("c1ccccc1")                   # doctest: +SKIP
    """
    t0 = time.time()
    log.info("autosolve.start", smiles=smiles, llm=llm, az_model=az_model)

    if clear_cache:
        clear_cache_for_molecule(smiles)

    result = _run_retrosynthesis(
        smiles,
        llm=llm,
        az_model=az_model,
        stability_flag=stability_flag,
        hallucination_check=hallucination_check,
        use_protecting_group_feature=use_protecting_group_feature,
    )

    elapsed = time.time() - t0
    log.info("autosolve.done", smiles=smiles, elapsed_s=round(elapsed, 2))

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(result, indent=4))
        log.info("autosolve.saved", path=str(output_path))

    return result


async def autosolve_async(
    smiles_list: list[str],
    *,
    llm: str = "anthropic/claude-3-7-sonnet-20250219:adv",
    az_model: str = "Pistachio_100+",
    stability_flag: str = "True",
    hallucination_check: str = "True",
    use_protecting_group_feature: bool = False,
    clear_cache: bool = True,
    output_dir: str | Path | None = None,
    max_concurrent: int = 6,
) -> list[dict[str, Any]]:
    """Run retrosynthesis on multiple molecules concurrently.

    This is the async replacement for ``prod_dfs_runner_async.py``.
    Each molecule is dispatched to a :class:`ThreadPoolExecutor`
    (because the underlying pipeline is synchronous) and at most
    *max_concurrent* molecules run in parallel.

    Parameters
    ----------
    smiles_list : list[str]
        SMILES strings of the target molecules.
    llm : str
        LLM model identifier.
    az_model : str
        AiZynthFinder model variant.
    stability_flag : str
        ``"True"`` to enable stability checking.
    hallucination_check : str
        ``"True"`` to enable hallucination checking.
    use_protecting_group_feature : bool
        Enable protecting-group masking in the LLM prompt.
    clear_cache : bool
        Clear cached results for each molecule before running.
    output_dir : str or Path or None
        If given, each result is written as
        ``<output_dir>/<sanitised_smiles>.json``.
    max_concurrent : int
        Maximum number of molecules processed in parallel.

    Returns
    -------
    list[dict[str, Any]]
        One result dict per input SMILES, in the same order.
        Failed molecules produce ``{"error": "<msg>", "smiles": "<smi>"}``.

    Examples
    --------
    >>> import asyncio                                            # doctest: +SKIP
    >>> from deepretro.utils.autosolve import autosolve_async     # doctest: +SKIP
    >>> results = asyncio.run(autosolve_async(["c1ccccc1"]))      # doctest: +SKIP
    """
    log.info(
        "autosolve_async.start",
        count=len(smiles_list),
        max_concurrent=max_concurrent,
    )

    semaphore = asyncio.Semaphore(max_concurrent)
    loop = asyncio.get_running_loop()

    async def _process(smi: str) -> dict[str, Any]:
        async with semaphore:
            t0 = time.time()
            try:
                if clear_cache:
                    clear_cache_for_molecule(smi)

                with ThreadPoolExecutor(max_workers=1) as pool:
                    result = await loop.run_in_executor(
                        pool,
                        functools.partial(
                            _run_retrosynthesis,
                            smi,
                            llm=llm,
                            az_model=az_model,
                            stability_flag=stability_flag,
                            hallucination_check=hallucination_check,
                            use_protecting_group_feature=use_protecting_group_feature,
                        ),
                    )

                if output_dir is not None:
                    out = Path(output_dir)
                    out.mkdir(parents=True, exist_ok=True)
                    fname = _safe_filename(smi) + ".json"
                    (out / fname).write_text(json.dumps(result, indent=4))

                elapsed = time.time() - t0
                log.info(
                    "autosolve_async.mol_done",
                    smiles=smi,
                    elapsed_s=round(elapsed, 2),
                )
                return result

            except Exception as exc:
                elapsed = time.time() - t0
                log.error(
                    "autosolve_async.mol_error",
                    smiles=smi,
                    error=str(exc),
                    elapsed_s=round(elapsed, 2),
                )
                return {"error": str(exc), "smiles": smi}

    tasks = [asyncio.create_task(_process(smi)) for smi in smiles_list]
    results = await asyncio.gather(*tasks)

    log.info("autosolve_async.done", total=len(results))
    return list(results)
