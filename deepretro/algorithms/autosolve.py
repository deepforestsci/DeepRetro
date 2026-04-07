"""Batch and single-molecule retrosynthesis runners.

Wraps the full DeepRetro pipeline (AiZynthFinder -> LLM fallback ->
metadata enrichment) into callable functions for headless execution
without the GUI or server.

The pipeline flow for each molecule is:

1. **AiZynthFinder** attempts a template-based retrosynthesis.
2. If AZ cannot find a route, the **LLM** (e.g. Claude) proposes
   retrosynthetic steps with chain-of-thought reasoning.
3. Proposed pathways are validated: SMILES validity, optional
   **stability checks** (ring strain, reactive intermediates),
   and optional **hallucination detection** (heuristic or ML-based).
4. The recursive pipeline repeats for each reactant until all
   leaf nodes are purchasable (in-stock) chemicals.
5. The result tree is flattened into a step-by-step synthesis plan
   with dependency ordering, confidence estimates, and metadata.

Public API
----------
* :func:`autosolve` — synchronous, one molecule at a time.
* :func:`autosolve_async` — concurrent execution over a list of SMILES
  using ``asyncio`` + a thread-pool (the underlying pipeline is
  synchronous, so each molecule gets its own worker thread).

Examples
--------
>>> from deepretro.algorithms.autosolve import autosolve  # doctest: +SKIP
>>> result = autosolve("CC(=O)Oc1ccccc1C(=O)O")           # doctest: +SKIP
>>> result["steps"][0]["reactants"]                         # doctest: +SKIP
[{'smiles': 'CC(=O)O'}, {'smiles': 'Oc1ccccc1C(=O)O'}]
"""

from __future__ import annotations

import asyncio
import functools
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from deepretro.models.hallucination_helpers import resolve_hallucination_args


def get_pipeline():
    """Lazy-import the retrosynthesis pipeline and return a runner callable.

    Imports are deferred so that ``import deepretro.algorithms.autosolve``
    succeeds even when heavy dependencies (``aizynthfinder``, LLM client
    libraries) are not installed — useful for CI, docs builds, and
    lightweight downstream packages that only need the data utilities.

    Returns
    -------
    callable
        A function ``run_retrosynthesis(smiles, **kwargs) -> dict`` that
        executes one full retrosynthesis pass and returns the formatted
        result dictionary.
    """
    import os
    import time

    import structlog
    from deepretro.migration.rec_prithvi import rec_run_prithvi
    from deepretro.migration.job_context import logger as context_logger
    from deepretro.migration.parse import format_output

    def run_retrosynthesis(smiles, **kwargs):
        """Execute one retrosynthesis and return the formatted result.

        Parameters
        ----------
        smiles : str
            Target molecule SMILES.
        **kwargs
            Forwarded to :func:`rec_run_prithvi` — includes ``llm``,
            ``az_model``, ``stability_flag``, ``hallucination_check``,
            ``use_protecting_group_feature``, ``hallucination_checker_fn``.

        Returns
        -------
        dict
            Formatted result with ``steps`` and ``dependencies`` keys.
        """
        job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
        log = structlog.get_logger().bind(job_id=job_id)
        token = context_logger.set(log)
        try:
            result_dict, _ = rec_run_prithvi(smiles, job_id=job_id, **kwargs)
            return format_output(result_dict)
        finally:
            context_logger.reset(token)

    return run_retrosynthesis


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

    This is the main entry point for programmatic use of DeepRetro.
    It is a pure function: returns the result dict without writing
    anything to disk.  The caller decides what to do with the output
    (save, display, post-process, etc.).

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.
    llm : str
        LLM model identifier passed to the pipeline.  Use the format
        ``"provider/model:adv:think"`` — the ``:adv`` suffix selects
        the advanced prompt template and ``:think`` enables extended
        thinking for supported models (e.g. Claude Opus 4.6).
    az_model : str
        AiZynthFinder model variant (e.g. ``"USPTO"``,
        ``"Pistachio_100+"``).  Must match a subdirectory under
        ``AZ_MODELS_PATH`` or default to the fallback config.
    stability_check : bool
        When ``True``, proposed reactants are checked for chemical
        stability (ring strain, carbocations, carbenes, etc.) and
        unstable pathways are discarded before returning results.
    hallucination_mode : str
        Which hallucination checker to use inside the pipeline's
        retry loop.  One of:

        * ``"heuristic"`` — rule-based checker that compares product
          and reactant atom counts, molecular weight, and fingerprint
          similarity (default).
        * ``"ml"`` — use the ML ``HallucinationClassifier`` supplied
          via *hallucination_classifier*.  Requires a fitted model.
        * ``"none"`` — skip hallucination checking entirely.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when ``hallucination_mode="ml"``.  Pass a fitted
        :class:`~deepretro.models.HallucinationClassifier` instance
        or a path to a saved model directory (str / Path).
    use_protecting_group_feature : bool
        When ``True``, the LLM prompt includes protecting-group context
        derived from automatic detection of common protecting groups
        in the target molecule.

    Returns
    -------
    dict[str, Any]
        Retrosynthesis result dictionary with keys:

        * ``"steps"`` — list of reaction steps, each containing
          ``reactants``, ``products``, ``conditions``, ``reagents``,
          and ``reactionmetrics`` (including ``confidenceestimate``).
        * ``"dependencies"`` — adjacency mapping from each step to
          its prerequisite steps.

    Raises
    ------
    ValueError
        If *hallucination_mode* is invalid, or ``"ml"`` mode is
        requested without a valid *hallucination_classifier*.

    Examples
    --------
    >>> from deepretro.algorithms.autosolve import autosolve  # doctest: +SKIP
    >>> result = autosolve("c1ccccc1")                        # doctest: +SKIP
    >>> result = autosolve(
    ...     "c1ccccc1",
    ...     hallucination_mode="ml",
    ...     hallucination_classifier="model_out/",
    ... )                                                      # doctest: +SKIP
    """
    run_retrosynthesis = get_pipeline()
    h_check, h_checker_fn = resolve_hallucination_args(
        hallucination_mode, hallucination_classifier,
    )

    return run_retrosynthesis(
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

    Each molecule is dispatched to a :class:`ThreadPoolExecutor`
    (because the underlying pipeline is synchronous) and at most
    *max_concurrent* molecules run in parallel.  Returns results
    in the same order as the input list.

    This is a pure function: returns a list of result dicts without
    writing anything to disk.

    Parameters
    ----------
    smiles_list : list[str]
        SMILES strings of the target molecules.
    llm : str
        LLM model identifier.  See :func:`autosolve` for format.
    az_model : str
        AiZynthFinder model variant.
    stability_check : bool
        Enable stability checking on proposed reactants.
    hallucination_mode : str
        ``"heuristic"``, ``"ml"``, or ``"none"``.
        See :func:`autosolve` for details.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when ``hallucination_mode="ml"``.
        See :func:`autosolve` for details.
    use_protecting_group_feature : bool
        Enable protecting-group masking in the LLM prompt.
    max_concurrent : int
        Maximum number of molecules processed in parallel.  Each
        molecule gets its own thread; this caps the concurrency to
        avoid overwhelming API rate limits.

    Returns
    -------
    list[dict[str, Any]]
        One result dict per input SMILES, in the same order.
        Failed molecules produce ``{"error": "<msg>", "smiles": "<smi>"}``
        instead of raising, so one failure does not abort the batch.

    Examples
    --------
    >>> import asyncio                                                 # doctest: +SKIP
    >>> from deepretro.algorithms.autosolve import autosolve_async     # doctest: +SKIP
    >>> results = asyncio.run(autosolve_async(["c1ccccc1", "CCO"]))    # doctest: +SKIP
    >>> len(results)                                                    # doctest: +SKIP
    2
    """
    run_retrosynthesis = get_pipeline()
    h_check, h_checker_fn = resolve_hallucination_args(
        hallucination_mode, hallucination_classifier,
    )

    stability_flag = str(stability_check)
    semaphore = asyncio.Semaphore(max_concurrent)
    loop = asyncio.get_running_loop()

    async def process_molecule(smi: str) -> dict[str, Any]:
        async with semaphore:
            try:
                with ThreadPoolExecutor(max_workers=1) as pool:
                    return await loop.run_in_executor(
                        pool,
                        functools.partial(
                            run_retrosynthesis,
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

    tasks = [asyncio.create_task(process_molecule(smi)) for smi in smiles_list]
    return list(await asyncio.gather(*tasks))
