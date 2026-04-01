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
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from deepretro.models.hallucination_helpers import resolve_hallucination_args


def _get_pipeline():
    """Lazy-import the retrosynthesis pipeline.

    Deferred so that ``import deepretro.utils.autosolve`` succeeds even
    when ``src`` / ``aizynthfinder`` are not installed (e.g. in CI).
    """
    import os
    import time

    import structlog
    from src.rec_prithvi import rec_run_prithvi
    from src.utils.job_context import logger as context_logger
    from src.utils.parse import format_output

    def _run(smiles, **kwargs):
        job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
        log = structlog.get_logger().bind(job_id=job_id)
        token = context_logger.set(log)
        try:
            result_dict, _ = rec_run_prithvi(smiles, job_id=job_id, **kwargs)
            return format_output(result_dict)
        finally:
            context_logger.reset(token)

    return _run


def autosolve(
    smiles: str,
    *,
    llm: str = "anthropic/claude-3-7-sonnet-20250219:adv",
    az_model: str = "Pistachio_100+",
    stability_check: bool = True,
    hallucination_mode: str = "heuristic",
    hallucination_classifier: Any = None,
    use_protecting_group_feature: bool = False,
) -> dict[str, Any]:
    """Run retrosynthesis on a single molecule.

    Pure function: returns the result dict without writing anything to
    disk.  The caller decides what to do with the output (save, display,
    post-process, etc.).

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.
    llm : str
        LLM model identifier passed to the pipeline.
    az_model : str
        AiZynthFinder model variant (e.g. ``"USPTO"``,
        ``"Pistachio_100+"``).
    stability_check : bool
        Enable stability checking on proposed reactants.
    hallucination_mode : str
        Which hallucination checker to use inside the pipeline's retry
        loop.  One of:

        * ``"heuristic"`` — rule-based checker (default, existing
          behaviour).
        * ``"ml"`` — use the ML ``HallucinationClassifier`` supplied
          via *hallucination_classifier*.
        * ``"none"`` — skip hallucination checking entirely.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when ``hallucination_mode="ml"``.  Pass a fitted
        :class:`~deepretro.models.HallucinationClassifier` instance
        or a path to a saved model directory.
    use_protecting_group_feature : bool
        Enable protecting-group masking in the LLM prompt.

    Returns
    -------
    dict[str, Any]
        Retrosynthesis result dictionary.

    Examples
    --------
    >>> from deepretro.utils.autosolve import autosolve  # doctest: +SKIP
    >>> result = autosolve("c1ccccc1")                   # doctest: +SKIP
    >>> result = autosolve("c1ccccc1", hallucination_mode="ml",
    ...     hallucination_classifier="saved_model/")      # doctest: +SKIP
    """
    _run_retrosynthesis = _get_pipeline()
    h_check, h_checker_fn = resolve_hallucination_args(
        hallucination_mode, hallucination_classifier,
    )

    return _run_retrosynthesis(
        smiles,
        llm=llm,
        az_model=az_model,
        stability_flag=str(stability_check),
        hallucination_check=h_check,
        use_protecting_group_feature=use_protecting_group_feature,
        hallucination_checker_fn=h_checker_fn,
    )


async def autosolve_async(
    smiles_list: list[str],
    *,
    llm: str = "anthropic/claude-3-7-sonnet-20250219:adv",
    az_model: str = "Pistachio_100+",
    stability_check: bool = True,
    hallucination_mode: str = "heuristic",
    hallucination_classifier: Any = None,
    use_protecting_group_feature: bool = False,
    max_concurrent: int = 6,
) -> list[dict[str, Any]]:
    """Run retrosynthesis on multiple molecules concurrently.

    Pure function: returns a list of result dicts without writing
    anything to disk.

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
    stability_check : bool
        Enable stability checking on proposed reactants.
    hallucination_mode : str
        ``"heuristic"``, ``"ml"``, or ``"none"``.
        See :func:`autosolve`.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when ``hallucination_mode="ml"``.
        See :func:`autosolve`.
    use_protecting_group_feature : bool
        Enable protecting-group masking in the LLM prompt.
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
    _run_retrosynthesis = _get_pipeline()
    h_check, h_checker_fn = resolve_hallucination_args(
        hallucination_mode, hallucination_classifier,
    )

    stability_flag = str(stability_check)
    semaphore = asyncio.Semaphore(max_concurrent)
    loop = asyncio.get_running_loop()

    async def _process(smi: str) -> dict[str, Any]:
        async with semaphore:
            try:
                with ThreadPoolExecutor(max_workers=1) as pool:
                    return await loop.run_in_executor(
                        pool,
                        functools.partial(
                            _run_retrosynthesis,
                            smi,
                            llm=llm,
                            az_model=az_model,
                            stability_flag=stability_flag,
                            hallucination_check=h_check,
                            use_protecting_group_feature=use_protecting_group_feature,
                            hallucination_checker_fn=h_checker_fn,
                        ),
                    )
            except Exception as exc:
                return {"error": str(exc), "smiles": smi}

    tasks = [asyncio.create_task(_process(smi)) for smi in smiles_list]
    return list(await asyncio.gather(*tasks))
