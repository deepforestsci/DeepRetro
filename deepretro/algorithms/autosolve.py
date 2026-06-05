"""Automatic single-molecule retrosynthesis solver."""

from __future__ import annotations

import os
import time
from collections.abc import Callable, Sequence
from typing import Any, Optional

import structlog

from deepretro.models.hallucination_helpers import (
    HallucinationPredictor,
    filter_with_checker,
    resolve_hallucination,
)
from deepretro.utils.az import run_az
from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.parse import format_output
from deepretro.metadata_types import (
    ConditionsRecommender,
    LiteratureRecommender,
    ReagentRecommender,
)
from deepretro.utils.typing import HallucinationChecker, ParseOutput, RouteNode
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
    hallucination_classifier : HallucinationPredictor or None, optional
        Loaded classifier object for ML hallucination filtering. Must expose
        ``predict_single(product_smiles, reactants_smiles)``.
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

    Raises
    ------
    ValueError
        If ``max_depth`` is negative, or if ``hallucination_mode`` is
        ``"ml"`` without a valid classifier.

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
        hallucination_classifier: HallucinationPredictor | None = None,
        max_depth: int = 50,
        enable_thinking: bool = True,
        max_output_tokens: int | None = None,
        metadata_model: str = "anthropic/claude-sonnet-4-6",
        az_runner: Callable[..., tuple[bool, list[Any]]] | None = None,
        llm_runner: Callable[..., tuple[list[Any], list[str], list[float]]]
        | None = None,
    ) -> None:
        """Construct an AutoSolver; see the class docstring for parameters."""
        if max_depth < 0:
            raise ValueError("max_depth must be non-negative")

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
        >>> route, solved = solver.solve("CC(=O)Oc1ccccc1C(=O)O")  # doctest: +SKIP
        >>> solved  # doctest: +SKIP
        False
        """
        smiles = _clean_smiles(smiles)
        if not smiles:
            return unsolved_leaf(smiles), False

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

        pathways, explanations, confidence = self.run_llm(smiles)
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
        >>> route, solved = solver.single_step("CC(=O)Oc1ccccc1C(=O)O")
        >>> solved
        False
        """
        smiles = _clean_smiles(smiles)
        if not smiles:
            return unsolved_leaf(smiles), False

        az_solved, az_routes = self._az_runner(smiles, self.az_model)
        if az_solved and az_routes:
            return dict(az_routes[0]), True

        pathways, explanations, confidence = self.run_llm(smiles)
        if not pathways:
            return unsolved_leaf(smiles), False

        for i, pathway in enumerate(pathways):
            reactants = _pathway_reactants(pathway)
            if not reactants:
                continue
            children = [unsolved_leaf(r) for r in reactants]
            selected_confidence = [confidence[i]] if i < len(confidence) else []
            return reaction_tree(smiles, children, selected_confidence), False

        return unsolved_leaf(smiles), False

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
        reagent_recommender: ReagentRecommender | None = None,
        conditions_recommender: ConditionsRecommender | None = None,
        literature_recommender: LiteratureRecommender | None = None,
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
        >>> output = solver.autosolve("CC(=O)Oc1ccccc1C(=O)O")  # doctest: +SKIP
        >>> output["solved"]  # doctest: +SKIP
        False
        """
        smiles = _clean_smiles(smiles)
        log = logger.bind(
            job_id=(f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}_{time.time_ns()}")
        )
        log.info("AutoSolver starting", molecule=smiles)

        route_tree, solved = self.solve(smiles)
        output = self.parse(route_tree, solved=solved)
        output = self.add_metadata(output)

        log.info("AutoSolver completed", molecule=smiles, solved=solved)
        return output

    def run_llm(
        self,
        molecule: str,
    ) -> tuple[list[Pathway], list[str], list[float]]:
        """Call the LLM retrosynthesis pipeline and filter hallucinations.

        Uses the injected ``llm_runner`` if provided at construction time,
        otherwise imports and calls ``deepretro.utils.llm.llm_pipeline``.
        Results are then passed through the hallucination checker (if any).

        Parameters
        ----------
        molecule : str
            Target molecule SMILES for the LLM pipeline.

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
