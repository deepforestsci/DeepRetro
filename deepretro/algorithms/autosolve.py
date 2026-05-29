"""Automatic single-molecule retrosynthesis solver."""

from __future__ import annotations

import os
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Optional

import structlog

from deepretro.models.hallucination_helpers import (
    HallucinationChecker,
    filter_with_checker,
    resolve_hallucination,
)
from deepretro.utils.az import run_az
from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.parse import format_output
from deepretro.utils.typing import ParseOutput, RouteNode
from deepretro.utils.utils_molecule import canonicalize

logger = structlog.get_logger(__name__)


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
        Hallucination filter strategy. The heuristic mode delegates to the
        package LLM pipeline; ML mode applies the provided classifier after
        the LLM pipeline returns candidates.
    hallucination_classifier : str, Path, object, or None, optional
        Saved classifier directory or classifier object for ML hallucination
        filtering.
    max_depth : int, optional
        Maximum retrosynthesis recursion depth.
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
        Replacement for the default LLM pipeline call. Accepts
        ``(molecule, **kwargs)`` and returns ``(pathways, explanations,
        confidence)``.

    Examples
    --------
    >>> solver = AutoSolver(hallucination_mode="none")
    >>> # Full pipeline (requires LLM and AZ services):
    >>> # result = solver.autosolve("CC(=O)Oc1ccccc1C(=O)O")
    >>> # Individual steps:
    >>> # route_tree, solved = solver.solve("CC(=O)Oc1ccccc1C(=O)O")
    >>> # parsed = solver.parse(route_tree, solved=solved)
    >>> # enriched = solver.add_metadata(parsed)
    """

    def __init__(
        self,
        llm: str = "anthropic/claude-opus-4-6:adv",
        az_model: str = "Pistachio_100+",
        stability_check: bool = True,
        hallucination_mode: str = "heuristic",
        hallucination_classifier: str | Path | Any | None = None,
        max_depth: int = 50,
        enable_thinking: bool = True,
        max_output_tokens: int | None = None,
        metadata_model: str = "claude-opus-4-20250514",
        az_runner: Callable[..., tuple[bool, list[Any]]] | None = None,
        llm_runner: Callable[..., tuple[list[Any], list[str], list[float]]]
        | None = None,
    ) -> None:
        self.llm = llm
        self.az_model = az_model
        self.stability_check = stability_check
        self.hallucination_mode = hallucination_mode.lower()
        self.hallucination_checker: Optional[HallucinationChecker] = (
            resolve_hallucination(
                self.hallucination_mode,
                hallucination_classifier,
            )
        )
        self.max_depth = max_depth
        self.enable_thinking = enable_thinking
        self.max_output_tokens = max_output_tokens
        self.metadata_model = metadata_model
        self._az_runner = az_runner if az_runner is not None else run_az
        self._llm_runner = llm_runner

    def solve(
        self,
        smiles: str,
        visited: set[str] | None = None,
        depth: int = 0,
    ) -> tuple[dict[str, Any], bool]:
        """Recursively solve one molecule via retrosynthesis.

        Attempts AiZynthFinder first. On failure, falls back to the LLM
        pipeline and recursively solves each proposed precursor.

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
        >>> route, solved = solver.solve("CC(=O)Oc1ccccc1C(=O)O")
        >>> solved
        False
        """
        canonical = canonicalize(smiles)
        if depth >= self.max_depth:
            logger.warning("Maximum recursion depth reached", molecule=smiles)
            return unsolved_leaf(smiles), False

        branch_visited = set() if visited is None else set(visited)
        if canonical in branch_visited:
            logger.warning("Cycle detected in retrosynthesis tree", molecule=smiles)
            return unsolved_leaf(smiles), False
        branch_visited.add(canonical)

        az_solved, az_routes = self._az_runner(smiles, self.az_model)
        if az_solved and az_routes:
            logger.info("AiZynthFinder solved molecule", molecule=smiles)
            return dict(az_routes[0]), True

        pathways, explanations, confidence = self._run_llm_fallback(smiles)
        if not pathways:
            logger.info("LLM fallback returned no usable pathways", molecule=smiles)
            return unsolved_leaf(smiles), False

        first_attempted_children: list[dict[str, Any]] | None = None
        for pathway in pathways:
            candidate_children, candidate_solved = self._solve_pathway(
                pathway,
                branch_visited,
                depth,
            )
            if first_attempted_children is None:
                first_attempted_children = candidate_children
            if candidate_solved:
                return reaction_tree(smiles, candidate_children, confidence), True

        return reaction_tree(
            smiles,
            [] if first_attempted_children is None else first_attempted_children,
            confidence,
        ), False

    def _run_llm_fallback(
        self,
        molecule: str,
    ) -> tuple[list[Pathway], list[str], list[float]]:
        """Call the LLM pipeline and apply optional ML hallucination filtering.

        Uses the injected ``llm_runner`` if provided, otherwise imports and
        calls ``deepretro.utils.llm.llm_pipeline`` directly.

        Parameters
        ----------
        molecule : str
            Target molecule SMILES for the LLM pipeline.

        Returns
        -------
        tuple[list[Pathway], list[str], list[float]]
            Filtered pathways, explanations, and confidence scores.
        """
        if self._llm_runner is not None:
            pathways, explanations, confidence = self._llm_runner(
                molecule,
                model=self.llm,
                stability_check=self.stability_check,
                hallucination_check=self.hallucination_mode == "heuristic",
                enable_thinking=self.enable_thinking,
                max_output_tokens=self.max_output_tokens,
            )
        else:
            from deepretro.utils.llm import llm_pipeline

            pathways, explanations, confidence = llm_pipeline(
                molecule=molecule,
                model=self.llm,
                stability_check=self.stability_check,
                hallucination_check=self.hallucination_mode == "heuristic",
                enable_thinking=self.enable_thinking,
                max_output_tokens=self.max_output_tokens,
            )
        return filter_with_checker(
            molecule,
            pathways,
            explanations,
            confidence,
            self.hallucination_checker,
        )

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
        reactants = [pathway] if isinstance(pathway, str) else list(pathway)
        children: list[dict[str, Any]] = []
        all_solved = True
        for reactant in reactants:
            child, solved = self.solve(reactant, visited, depth + 1)
            children.append(child)
            all_solved = all_solved and solved
        return children, all_solved

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
        >>> route, solved = solver.single_step("CC(=O)Oc1ccccc1C(=O)O")
        >>> solved
        False
        """
        az_solved, az_routes = self._az_runner(smiles, self.az_model)
        if az_solved and az_routes:
            return dict(az_routes[0]), True

        pathways, explanations, confidence = self._run_llm_fallback(smiles)
        if not pathways:
            return unsolved_leaf(smiles), False

        first_pathway = pathways[0]
        reactants = (
            [first_pathway] if isinstance(first_pathway, str) else list(first_pathway)
        )
        children = [unsolved_leaf(r) for r in reactants]
        return reaction_tree(smiles, children, confidence), False

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
        return output

    def add_metadata(
        self,
        parsed_output: ParseOutput,
        *,
        reagent_recommender: Callable[..., Any] | None = None,
        conditions_recommender: Callable[..., Any] | None = None,
        literature_recommender: Callable[..., Any] | None = None,
    ) -> ParseOutput:
        """Enrich parsed steps with reagent, condition, and literature metadata.

        Iterates each step, builds a ``reactants>>product`` reaction SMILES,
        and calls :func:`deepretro.metadata.recommend_reaction_metadata`.
        Injectable recommender callables are forwarded, allowing tests to
        supply lightweight doubles.

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
        >>> from deepretro.algorithms.autosolve import reaction_tree, unsolved_leaf
        >>> solver = AutoSolver(hallucination_mode="none")
        >>> tree = reaction_tree("CCO", [unsolved_leaf("C"), unsolved_leaf("O")], [0.8])
        >>> parsed = solver.parse(tree, solved=True)
        >>> # enriched = solver.add_metadata(parsed)  # doctest: +SKIP
        """
        from deepretro.metadata import recommend_reaction_metadata

        for step in parsed_output.get("steps", []):
            reactant_smiles = ".".join(r["smiles"] for r in step.get("reactants", []))
            product_smiles = (
                step["products"][0]["smiles"] if step.get("products") else ""
            )
            if not reactant_smiles or not product_smiles:
                continue

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
                    step=step["step"],
                    status=status,
                )
                continue

            step["reagents"] = result.get("reagents", [])
            step["conditions"] = result.get("conditions", {})
            step["reactionmetrics"][0]["closestliterature"] = result.get(
                "literature", ""
            )

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
        >>> solver = AutoSolver(hallucination_mode="none")
        >>> # result = solver.autosolve("CC(=O)Oc1ccccc1C(=O)O")  # doctest: +SKIP
        >>> # print(result["solved"], len(result["steps"]))
        """
        log = logger.bind(job_id=f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}")
        log.info("AutoSolver starting", molecule=smiles)

        route_tree, solved = self.solve(smiles)
        output = self.parse(route_tree, solved=solved)
        output = self.add_metadata(output)

        log.info("AutoSolver completed", molecule=smiles, solved=solved)
        return output


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
    >>> unsolved_leaf("CC(=O)Oc1ccccc1C(=O)O")["in_stock"]
    False
    """
    return {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": False,
        "children": [],
    }


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
