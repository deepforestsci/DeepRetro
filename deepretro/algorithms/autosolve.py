"""Automatic single-molecule retrosynthesis solver.

Ports the recursive AiZynthFinder + LLM retrosynthesis loop into the
``deepretro`` package. AZ is tried first at every node; when it cannot solve a
molecule the LLM (or, under ``solve_mode="single_step_agent"``, a tool-calling
agent) proposes precursors, which are recursively solved until AZ solves each
leaf. ``run_llm`` is the single owner of validity, stability, and hallucination
filtering, so the agent path can never bypass the safety checks.
"""

from __future__ import annotations

import os
import time
from collections.abc import Callable, Sequence
from typing import Any, Optional

import structlog

from deepretro.metadata_types import (
    ConditionsRecommender,
    LiteratureRecommender,
    ReagentRecommender,
)
from deepretro.models.hallucination_helpers import (
    filter_with_checker,
    resolve_hallucination,
)
from deepretro.utils.az import run_az
from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.parse import format_output
from deepretro.utils.typing import HallucinationChecker, ParseOutput, RouteNode
from deepretro.utils.utils_molecule import canonicalize, validity_check

logger = structlog.get_logger(__name__)

VALID_SOLVE_MODES = ("pipeline", "single_step_agent", "orchestrator")
VALID_TOOL_BACKENDS = ("structured", "sandbox")


class AutoSolver:
    """Solve a target molecule with AiZynthFinder and an LLM fallback.

    Parameters
    ----------
    llm : str, optional
        LiteLLM model identifier used by the fallback retrosynthesis pipeline.
    az_model : str, optional
        AiZynthFinder model variant.
    stability_check : bool, optional
        Whether LLM candidates should pass the heuristic stability filter.
    hallucination_mode : {"heuristic", "ml", "none"}, optional
        Hallucination filter strategy. ``heuristic`` and ``ml`` both resolve to
        a checker that runs inside :meth:`run_llm`; ``none`` disables filtering.
    hallucination_classifier : object or str or Path or None, optional
        For ``hallucination_mode="ml"``, either a loaded classifier exposing
        ``predict_probability``/``threshold``/``featurizer`` or a path to a
        saved model directory. Ignored for the other modes.
    solve_mode : {"pipeline", "single_step_agent", "orchestrator"}, optional
        ``pipeline`` (default) uses the non-agentic LLM pipeline. ``single_step_agent``
        lets a tool-calling agent propose and self-check precursors. ``orchestrator``
        is reserved and raises :class:`NotImplementedError` when invoked.
    tool_backend : {"structured", "sandbox"}, optional
        Which tools the agent may call. Ignored when ``solve_mode="pipeline"``.
    sandbox : Sandbox or None, optional
        Sandbox used by the agent's ``run_python`` tool. Injectable for testing.
    agent_az_tools : bool, optional
        When ``True`` (and ``solve_mode="single_step_agent"``), expose the
        ``az_suggest`` / ``az_route`` AiZynthFinder tools to the agent so it can
        consult AZ for template-grounded disconnections. Ignored otherwise.
    max_depth : int, optional
        Maximum retrosynthesis recursion depth. This is a safety guard against
        runaway recursion; reaching it is logged as a warning.
    stop_depth : int or None, optional
        User-defined recursion cutoff. When set, any molecule reached at
        ``depth >= stop_depth`` is returned as an ordinary unsolved leaf
        without being expanded (no AZ or LLM call) and without the
        ``max_depth`` warning. ``None`` (default) disables the cutoff, leaving
        ``max_depth`` as the only limit. Independent of ``max_depth``, which
        still guards against runaway recursion.
    enable_thinking : bool, optional
        Whether provider-supported reasoning controls should be enabled.
    max_output_tokens : int, optional
        Output-token override for LLM calls.
    metadata_model : str, optional
        Model identifier used for metadata recommendation steps.
    az_runner : callable, optional
        Replacement for the default AiZynthFinder call. Accepts
        ``(smiles, az_model)`` and returns ``(solved, routes)``.
    llm_runner : callable, optional
        Replacement for the default LLM pipeline call (``pipeline`` mode).
    agent_runner : callable, optional
        Replacement for the default tool-calling agent (``single_step_agent`` mode).
        Accepts ``(molecule, **kwargs)`` and returns ``(pathways, explanations,
        confidence)``.

    Raises
    ------
    ValueError
        If ``max_depth`` or ``stop_depth`` is negative,
        ``solve_mode``/``tool_backend`` is unknown, or ``hallucination_mode``
        is ``"ml"`` without a classifier.

    Examples
    --------
    >>> solver = AutoSolver(
    ...     az_runner=lambda s, m: (False, []),
    ...     llm_runner=lambda mol, **kw: ([], [], []),
    ...     hallucination_mode="none",
    ... )
    >>> solver.max_depth
    50
    """

    def __init__(
        self,
        llm: str = "anthropic/claude-sonnet-4-6",
        az_model: str = "Pistachio_100+",
        stability_check: bool = True,
        hallucination_mode: str = "heuristic",
        hallucination_classifier: Any | None = None,
        solve_mode: str = "pipeline",
        tool_backend: str = "structured",
        sandbox: Any | None = None,
        agent_az_tools: bool = False,
        max_depth: int = 50,
        stop_depth: int | None = None,
        enable_thinking: bool = True,
        max_output_tokens: int | None = None,
        metadata_model: str = "anthropic/claude-sonnet-4-6",
        az_runner: Callable[..., tuple[bool, list[Any]]] | None = None,
        llm_runner: Callable[..., tuple[list[Any], list[str], list[float]]]
        | None = None,
        agent_runner: Callable[..., tuple[list[Any], list[str], list[float]]]
        | None = None,
    ) -> None:
        """Construct an AutoSolver; see the class docstring for parameters."""
        if max_depth < 0:
            raise ValueError("max_depth must be non-negative")
        if stop_depth is not None and stop_depth < 0:
            raise ValueError("stop_depth must be non-negative or None")

        solve_mode = solve_mode.lower()
        if solve_mode not in VALID_SOLVE_MODES:
            raise ValueError(
                f"solve_mode must be one of {VALID_SOLVE_MODES} — got {solve_mode!r}"
            )
        tool_backend = tool_backend.lower()
        if tool_backend not in VALID_TOOL_BACKENDS:
            raise ValueError(
                f"tool_backend must be one of {VALID_TOOL_BACKENDS} "
                f"— got {tool_backend!r}"
            )

        self.llm = llm
        self.az_model = az_model
        self.stability_check = stability_check
        self.hallucination_mode = hallucination_mode.lower()
        self.hallucination_checker: Optional[HallucinationChecker] = (
            resolve_hallucination(self.hallucination_mode, hallucination_classifier)
        )
        self.solve_mode = solve_mode
        self.tool_backend = tool_backend
        self.sandbox = sandbox
        self.agent_az_tools = agent_az_tools
        self.max_depth = max_depth
        self.stop_depth = stop_depth
        self.enable_thinking = enable_thinking
        self.max_output_tokens = max_output_tokens
        self.metadata_model = metadata_model
        self._az_runner = az_runner if az_runner is not None else run_az
        self._llm_runner = llm_runner
        self._agent_runner = agent_runner
        # Accumulates agent refusal / no-answer events across a solve run.
        self._agent_events: list[dict[str, Any]] = []

    # ------------------------------------------------------------------
    # Core pipeline methods
    # ------------------------------------------------------------------

    def solve(
        self,
        smiles: str,
        visited: set[str] | None = None,
        depth: int = 0,
    ) -> tuple[dict[str, Any], bool]:
        """Recursively solve one molecule via retrosynthesis.

        Attempts AiZynthFinder first. On failure, falls back to the LLM
        pipeline (or agent) and recursively solves each proposed precursor
        until AZ solves every leaf.

        Parameters
        ----------
        smiles : str
            Molecule SMILES to solve.
        visited : set[str], optional
            Canonical SMILES already visited on the current branch.
        depth : int, optional
            Current recursion depth.

        Returns
        -------
        tuple[dict[str, Any], bool]
            Raw route tree and solved flag.

        Examples
        --------
        >>> solver = AutoSolver(
        ...     az_runner=lambda s, m: (False, []),
        ...     llm_runner=lambda mol, **kw: ([], [], []),
        ...     hallucination_mode="none",
        ... )
        >>> route, solved = solver.solve("CCO")  # doctest: +SKIP
        >>> solved  # doctest: +SKIP
        False
        """
        self._reject_orchestrator()
        smiles = _clean_smiles(smiles)
        if not smiles:
            return unsolved_leaf(smiles), False

        canonical = canonicalize(smiles)
        if self.stop_depth is not None and depth >= self.stop_depth:
            logger.debug(
                "Reached configured stop_depth; truncating branch",
                molecule=smiles,
                stop_depth=self.stop_depth,
            )
            return unsolved_leaf(smiles), False
        if depth >= self.max_depth:
            logger.warning("Maximum recursion depth reached", molecule=smiles)
            return unsolved_leaf(smiles), False

        branch_visited = set() if visited is None else set(visited)
        if canonical in branch_visited:
            logger.warning("Cycle detected in retrosynthesis tree", molecule=smiles)
            return unsolved_leaf(smiles), False
        branch_visited.add(canonical)

        az_solved, az_routes, az_errored = self._run_az(smiles)
        if az_errored:
            return unsolved_leaf(smiles), False
        if az_solved and az_routes:
            logger.info("AiZynthFinder solved molecule", molecule=smiles)
            route = dict(az_routes[0])
            _mark_az_generated(route)
            return route, True

        pathways, _explanations, confidence = self.run_llm(smiles)
        if not pathways:
            logger.info("LLM fallback returned no usable pathways", molecule=smiles)
            return unsolved_leaf(smiles), False

        first_attempted_children: list[dict[str, Any]] | None = None
        first_attempted_confidence: list[float] = []
        for i, pathway in enumerate(pathways):
            candidate_children, candidate_solved = self._solve_pathway(
                pathway,
                branch_visited,
                depth,
            )
            if not candidate_children:
                continue
            if first_attempted_children is None:
                first_attempted_children = candidate_children
                first_attempted_confidence = (
                    [confidence[i]] if i < len(confidence) else []
                )
            if candidate_solved:
                selected_confidence = [confidence[i]] if i < len(confidence) else []
                return (
                    reaction_tree(smiles, candidate_children, selected_confidence),
                    True,
                )

        if first_attempted_children is None:
            return unsolved_leaf(smiles), False

        return reaction_tree(
            smiles, first_attempted_children, first_attempted_confidence
        ), False

    def single_step(
        self,
        smiles: str,
    ) -> tuple[dict[str, Any], bool]:
        """Run a single AZ + LLM retrosynthesis pass without recursion.

        Unlike :meth:`solve`, this does not recurse into the proposed
        precursors. Each precursor is returned as an unsolved leaf.

        Parameters
        ----------
        smiles : str
            Target molecule SMILES.

        Returns
        -------
        tuple[dict[str, Any], bool]
            Route tree (one level deep) and solved flag.

        Examples
        --------
        >>> solver = AutoSolver(
        ...     az_runner=lambda s, m: (False, []),
        ...     llm_runner=lambda mol, **kw: ([], [], []),
        ...     hallucination_mode="none",
        ... )
        >>> route, solved = solver.single_step("CCO")
        >>> solved
        False
        """
        self._reject_orchestrator()
        smiles = _clean_smiles(smiles)
        if not smiles:
            return unsolved_leaf(smiles), False

        az_solved, az_routes, az_errored = self._run_az(smiles)
        if az_errored:
            return unsolved_leaf(smiles), False
        if az_solved and az_routes:
            return dict(az_routes[0]), True

        pathways, _explanations, confidence = self.run_llm(smiles)
        if not pathways:
            return unsolved_leaf(smiles), False

        for i, pathway in enumerate(pathways):
            reactants = _pathway_reactants(pathway)
            if not reactants:
                continue
            children = [unsolved_leaf(reactant) for reactant in reactants]
            selected_confidence = [confidence[i]] if i < len(confidence) else []
            return reaction_tree(smiles, children, selected_confidence), False

        return unsolved_leaf(smiles), False

    def solve_multiple(
        self,
        smiles: str,
        k: int = 3,
    ) -> list[tuple[dict[str, Any], bool]]:
        """Solve the top-``k`` candidate first-step pathways independently.

        When AZ solves the target directly, a single solved route is returned.
        Otherwise each of the first ``k`` LLM/agent pathways is solved with its
        own fresh visited set, so shared reactants do not cross-contaminate
        cycle detection. This produces the multiple routes the batch runner
        writes as ``pathway_<i>.json``.

        Parameters
        ----------
        smiles : str
            Target molecule SMILES.
        k : int, optional
            Maximum number of candidate routes to return.

        Returns
        -------
        list[tuple[dict[str, Any], bool]]
            Up to ``k`` ``(route_tree, solved)`` tuples (at least one).

        Examples
        --------
        >>> solver = AutoSolver(
        ...     az_runner=lambda s, m: (True, [{"type": "mol", "smiles": s,
        ...         "is_chemical": True, "in_stock": True}]),
        ...     hallucination_mode="none",
        ... )
        >>> [solved for _, solved in solver.solve_multiple("CCO")]
        [True]
        """
        self._reject_orchestrator()
        smiles = _clean_smiles(smiles)
        if not smiles:
            return [(unsolved_leaf(smiles), False)]

        az_solved, az_routes, az_errored = self._run_az(smiles)
        if az_errored:
            return [(unsolved_leaf(smiles), False)]
        if az_solved and az_routes:
            return [(dict(az_routes[0]), True)]

        pathways, _explanations, confidence = self.run_llm(smiles)
        if not pathways:
            return [(unsolved_leaf(smiles), False)]

        canonical = canonicalize(smiles)
        results: list[tuple[dict[str, Any], bool]] = []
        for i, pathway in enumerate(pathways[:k]):
            # Fresh visited set per route so a reactant used by one candidate
            # route is not treated as a cycle in another.
            candidate_children, candidate_solved = self._solve_pathway(
                pathway,
                {canonical},
                depth=0,
            )
            if not candidate_children:
                continue
            selected_confidence = [confidence[i]] if i < len(confidence) else []
            results.append(
                (
                    reaction_tree(smiles, candidate_children, selected_confidence),
                    candidate_solved,
                )
            )

        if not results:
            return [(unsolved_leaf(smiles), False)]
        return results

    def parse(self, route_tree: dict[str, Any], *, solved: bool) -> ParseOutput:
        """Format a route tree into viewer-ready steps and dependencies.

        Delegates to :func:`deepretro.utils.parse.format_output` and appends
        the ``solved`` flag.

        Parameters
        ----------
        route_tree : dict[str, Any]
            Route tree from :meth:`solve` or :meth:`single_step`.
        solved : bool
            Whether the route was fully solved.

        Returns
        -------
        dict[str, Any]
            Parsed output with ``steps``, ``dependencies``, and ``solved`` keys.

        Examples
        --------
        >>> from deepretro.algorithms.autosolve import unsolved_leaf
        >>> solver = AutoSolver(hallucination_mode="none")
        >>> result = solver.parse(unsolved_leaf("CCO"), solved=False)
        >>> result["solved"]
        False
        """
        output = format_output(route_tree)
        output["solved"] = solved
        summary = summarize_az(route_tree)
        output["az_solved"] = summary["az_solved"]
        output["az_summary"] = summary
        return output

    def add_metadata(
        self,
        parsed_output: ParseOutput,
        *,
        reagent_recommender: ReagentRecommender | None = None,
        conditions_recommender: ConditionsRecommender | None = None,
        literature_recommender: LiteratureRecommender | None = None,
    ) -> ParseOutput:
        """Enrich parsed steps with reagent, condition, and literature metadata.

        Iterates each step, builds a ``reactants>>product`` reaction SMILES,
        and calls :func:`deepretro.metadata.recommend_reaction_metadata`.

        Parameters
        ----------
        parsed_output : dict[str, Any]
            Output from :meth:`parse`.
        reagent_recommender : callable, optional
            Override the default reagent prediction agent.
        conditions_recommender : callable, optional
            Override the default conditions prediction agent.
        literature_recommender : callable, optional
            Override the default literature prediction agent.

        Returns
        -------
        dict[str, Any]
            The same dict with each step enriched in-place.

        Examples
        --------
        >>> solver = AutoSolver(hallucination_mode="none")
        >>> solver.add_metadata({"steps": [], "dependencies": {}, "solved": True})
        {'steps': [], 'dependencies': {}, 'solved': True}
        """
        for step in parsed_output.get("steps", []):
            reactant_smiles = ".".join(r["smiles"] for r in step.get("reactants", []))
            product_smiles = (
                step["products"][0]["smiles"] if step.get("products") else ""
            )
            if not reactant_smiles or not product_smiles:
                continue

            from deepretro.metadata import recommend_reaction_metadata

            reaction_smiles = f"{reactant_smiles}>>{product_smiles}"
            status, result = recommend_reaction_metadata(
                reaction_smiles,
                model=self.metadata_model,
                reagent_recommender=reagent_recommender,
                conditions_recommender=conditions_recommender,
                literature_recommender=literature_recommender,
                cache=None,
            )
            if status != 200:
                logger.info(
                    "metadata recommendation failed",
                    step=step.get("step"),
                    status=status,
                )
                continue

            step["reagents"] = result.get("reagents", [])
            step["conditions"] = result.get("conditions", {})
            reactionmetrics = step.setdefault("reactionmetrics", [])
            if not reactionmetrics:
                reactionmetrics.append({"scalabilityindex": 0, "closestliterature": ""})
            # ``literature`` is the recommender payload (typically a
            # ``{"doi": ...}`` mapping); ``closestliterature`` is a string DOI.
            literature = result.get("literature", "")
            if isinstance(literature, dict):
                literature = literature.get("doi", "")
            reactionmetrics[0]["closestliterature"] = literature

        return parsed_output

    def autosolve(self, smiles: str) -> ParseOutput:
        """Run the full retrosynthesis pipeline: solve, parse, and enrich.

        Parameters
        ----------
        smiles : str
            Target molecule SMILES.

        Returns
        -------
        dict[str, Any]
            Fully enriched output with steps, dependencies, metadata, and
            solved status.

        Examples
        --------
        >>> solver = AutoSolver(
        ...     az_runner=lambda s, m: (False, []),
        ...     llm_runner=lambda mol, **kw: ([], [], []),
        ...     hallucination_mode="none",
        ... )
        >>> output = solver.autosolve("CCO")  # doctest: +SKIP
        >>> output["solved"]  # doctest: +SKIP
        False
        """
        smiles = _clean_smiles(smiles)
        log = logger.bind(
            job_id=(f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}_{time.time_ns()}")
        )
        log.info("AutoSolver starting", molecule=smiles)

        self._agent_events = []
        route_tree, solved = self.solve(smiles)
        output = self.parse(route_tree, solved=solved)
        output = self.add_metadata(output)
        if self._agent_events:
            output["agent_events"] = list(self._agent_events)
            output["agent_refusal"] = any(
                event.get("kind") == "refusal" for event in self._agent_events
            )

        log.info("AutoSolver completed", molecule=smiles, solved=solved)
        return output

    def run_llm(
        self,
        molecule: str,
    ) -> tuple[list[Pathway], list[str], list[float]]:
        """Obtain candidate pathways and apply all safety filters.

        This is the single owner of validity, stability, and hallucination
        filtering. In ``pipeline`` mode the LLM pipeline already applies
        validity and stability, so only hallucination is applied here. In
        ``single_step_agent`` mode the raw agent output is passed through the
        full validity → stability → hallucination filter chain, so the agent
        can never bypass the safety net.

        Parameters
        ----------
        molecule : str
            Target molecule SMILES for the LLM/agent call.

        Returns
        -------
        tuple[list[Pathway], list[str], list[float]]
            Filtered pathways, explanations, and confidence scores.

        Examples
        --------
        >>> solver = AutoSolver(
        ...     llm_runner=lambda mol, **kw: ([["CCO"]], ["reduction"], [0.8]),
        ...     hallucination_mode="none",
        ... )
        >>> pathways, explanations, confidence = solver.run_llm("CC=O")
        >>> pathways
        [['CCO']]
        """
        if self.solve_mode == "single_step_agent":
            pathways, explanations, confidence = self._run_agent(molecule)
            pathways, explanations, confidence = self._apply_safety_filters(
                molecule, pathways, explanations, confidence
            )
        else:
            pathways, explanations, confidence = self._run_pipeline(molecule)

        return filter_with_checker(
            molecule,
            pathways,
            explanations,
            confidence,
            self.hallucination_checker,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _reject_orchestrator(self) -> None:
        """Raise for the not-yet-implemented top-level orchestrator mode."""
        if self.solve_mode == "orchestrator":
            raise NotImplementedError(
                "solve_mode='orchestrator' is not implemented yet; use "
                "'pipeline' or 'single_step_agent'"
            )

    def _run_az(self, smiles: str) -> tuple[bool, list[Any], bool]:
        """Call the AZ runner, converting any failure into an error flag.

        Returns ``(solved, routes, errored)``. AiZynthFinder is an external,
        occasionally-misconfigured tool (missing model files raise
        ``FileNotFoundError``/``ImportError``); a broad guard keeps one node's
        failure from crashing the whole solve.
        """
        try:
            solved, routes = self._az_runner(smiles, self.az_model)
            return solved, list(routes), False
        except Exception as exc:  # ponytail: AZ is external; never crash a solve
            logger.error("AZ runner failed", molecule=smiles, error=str(exc))
            return False, [], True

    def _run_pipeline(
        self, molecule: str
    ) -> tuple[list[Pathway], list[str], list[float]]:
        """Get raw pathways from the non-agentic LLM pipeline."""
        runner = self._llm_runner
        if runner is None:
            from deepretro.utils.llm import llm_pipeline

            runner = llm_pipeline
        return runner(
            molecule,
            model=self.llm,
            stability_check=self.stability_check,
            hallucination_check=False,
            enable_thinking=self.enable_thinking,
            max_output_tokens=self.max_output_tokens,
        )

    def _run_agent(self, molecule: str) -> tuple[list[Pathway], list[str], list[float]]:
        """Get raw pathways from the tool-calling agent."""
        runner = self._agent_runner
        kwargs: dict[str, Any] = dict(
            model=self.llm,
            tool_backend=self.tool_backend,
            sandbox=self.sandbox,
            hallucination_checker=self.hallucination_checker,
            enable_thinking=self.enable_thinking,
            max_output_tokens=self.max_output_tokens,
        )
        if runner is None:
            from deepretro.agents.loop import agentic_single_step

            runner = agentic_single_step
            # Only the built-in runner accepts these; a custom runner may not.
            kwargs["event_sink"] = self._agent_events
            kwargs["az_tools"] = self.agent_az_tools
            kwargs["az_model"] = self.az_model
        return runner(molecule, **kwargs)

    def _apply_safety_filters(
        self,
        molecule: str,
        pathways: list[Pathway],
        explanations: list[str],
        confidence: list[float],
    ) -> tuple[list[Pathway], list[str], list[float]]:
        """Apply validity then optional stability filtering, keeping alignment."""
        pathways, explanations, confidence = validity_check(
            molecule, pathways, explanations, confidence
        )
        if self.stability_check:
            from deepretro.utils.llm import (
                filter_metadata_by_pathways,
                filter_stable_pathways,
            )

            stable = filter_stable_pathways(pathways)
            pathways, explanations, confidence = filter_metadata_by_pathways(
                stable, pathways, explanations, confidence
            )
        return pathways, explanations, confidence

    def _solve_pathway(
        self,
        pathway: Sequence[str] | str,
        visited: set[str],
        depth: int,
    ) -> tuple[list[dict[str, Any]], bool]:
        """Solve each reactant in a candidate precursor pathway.

        Parameters
        ----------
        pathway : Sequence[str] or str
            One or more reactant SMILES from a candidate retrosynthesis step.
        visited : set[str]
            Canonical SMILES already visited on the current branch.
        depth : int
            Current recursion depth.

        Returns
        -------
        tuple[list[dict[str, Any]], bool]
            Solved child route nodes and whether all reactants were solved.
        """
        reactants = _pathway_reactants(pathway)
        if not reactants:
            return [], False
        children: list[dict[str, Any]] = []
        all_solved = True
        for reactant in reactants:
            child, solved = self.solve(reactant, visited, depth + 1)
            children.append(child)
            all_solved = all_solved and solved
        return children, all_solved


# ----------------------------------------------------------------------
# Module-level route tree builders
# ----------------------------------------------------------------------


def unsolved_leaf(smiles: str) -> dict[str, Any]:
    """Build a terminal route node for an unresolved molecule.

    Parameters
    ----------
    smiles : str
        SMILES of the unresolved molecule.

    Returns
    -------
    dict[str, Any]
        Route node with ``in_stock=False`` and empty children.

    Examples
    --------
    >>> unsolved_leaf("CCO")["in_stock"]
    False
    """
    return {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": False,
        "children": [],
    }


def _mark_az_generated(node: dict[str, Any]) -> None:
    """Recursively mark every node of an AiZynthFinder-produced subtree.

    Sets ``az_generated=True`` on each node so that, in a mixed route tree
    (LLM proposals with AZ-solved sub-branches), the nodes AZ generated can be
    told apart from the LLM-proposed ones.

    Examples
    --------
    >>> node = {"type": "mol", "smiles": "CCO", "children": []}
    >>> _mark_az_generated(node)
    >>> node["az_generated"]
    True
    """
    if not isinstance(node, dict):
        return
    node["az_generated"] = True
    for child in node.get("children") or []:
        _mark_az_generated(child)


def _collect_leaf_nodes(node: dict[str, Any], acc: list[dict[str, Any]]) -> None:
    """Collect every leaf ``mol`` node (a ``mol`` node with no children)."""
    children = node.get("children") or []
    if node.get("type") == "mol" and not children:
        acc.append(node)
        return
    for child in children:
        _collect_leaf_nodes(child, acc)


def summarize_az(route_tree: dict[str, Any]) -> dict[str, Any]:
    """Summarise AiZynthFinder's contribution to a route tree.

    Parameters
    ----------
    route_tree : dict[str, Any]
        Raw route tree from :meth:`AutoSolver.solve`, with AZ-generated
        subtrees marked by :func:`_mark_az_generated`.

    Returns
    -------
    dict[str, Any]
        ``az_solved`` (AZ generated at least one leaf), ``az_solved_all``
        (AZ generated *every* leaf, i.e. it closed all terminal branches),
        and leaf counts (total, AZ-generated, in-stock, unsolved).

    Examples
    --------
    >>> leaf = {"type": "mol", "smiles": "CCO", "in_stock": True,
    ...         "az_generated": True, "children": []}
    >>> summarize_az(leaf)["az_solved_all"]
    True
    """
    leaves: list[dict[str, Any]] = []
    _collect_leaf_nodes(route_tree, leaves)
    total = len(leaves)
    az_generated = sum(1 for leaf in leaves if leaf.get("az_generated"))
    in_stock = sum(1 for leaf in leaves if leaf.get("in_stock"))
    return {
        "az_solved": az_generated > 0,
        "az_solved_all": total > 0 and az_generated == total,
        "leaves_total": total,
        "leaves_az_generated": az_generated,
        "leaves_in_stock": in_stock,
        "leaves_unsolved": total - in_stock,
    }


def _clean_smiles(smiles: str) -> str:
    """Validate and strip surrounding whitespace from a SMILES string.

    Parameters
    ----------
    smiles : str
        Raw SMILES input.

    Returns
    -------
    str
        The stripped SMILES string.

    Raises
    ------
    TypeError
        If *smiles* is not a string.

    Examples
    --------
    >>> _clean_smiles("  CCO  ")
    'CCO'
    """
    if not isinstance(smiles, str):
        raise TypeError("smiles must be a string")
    return smiles.strip()


def _pathway_reactants(pathway: Sequence[str] | str) -> list[str]:
    """Normalize a pathway into a list of non-empty reactant SMILES.

    Parameters
    ----------
    pathway : Sequence[str] or str
        A single reactant string or a sequence of reactant SMILES.

    Returns
    -------
    list[str]
        Stripped, non-empty reactant SMILES.

    Raises
    ------
    TypeError
        If a sequence element is not a string.

    Examples
    --------
    >>> _pathway_reactants(["CCO", " ", "O"])
    ['CCO', 'O']
    >>> _pathway_reactants("CC=O")
    ['CC=O']
    """
    if isinstance(pathway, str):
        return [pathway.strip()] if pathway.strip() else []

    reactants: list[str] = []
    for reactant in pathway:
        if not isinstance(reactant, str):
            raise TypeError("pathway reactants must be strings")
        reactant = reactant.strip()
        if reactant:
            reactants.append(reactant)
    return reactants


def reaction_tree(
    molecule: str,
    children: Sequence[RouteNode],
    confidence: Sequence[float],
) -> dict[str, Any]:
    """Build the route-tree shape consumed by ``RetrosynthesisRouteParser``.

    Parameters
    ----------
    molecule : str
        Parent molecule SMILES.
    children : Sequence[RouteNode]
        Solved child route nodes (reactants).
    confidence : Sequence[float]
        Confidence scores from the LLM pipeline; the first value is used as
        ``policy_probability``.

    Returns
    -------
    dict[str, Any]
        Nested route tree with a single reaction child node.

    Examples
    --------
    >>> tree = reaction_tree("CCO", [unsolved_leaf("C"), unsolved_leaf("O")], [0.8])
    >>> tree["smiles"]
    'CCO'
    >>> tree["children"][0]["metadata"]["policy_probability"]
    0.8
    """
    policy_probability = float(confidence[0]) if confidence else 0.0
    return {
        "type": "mol",
        "smiles": molecule,
        "is_chemical": True,
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "is_reaction": True,
                "metadata": {"policy_probability": policy_probability},
                "children": list(children),
            }
        ],
    }
