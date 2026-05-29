"""Automatic single-molecule retrosynthesis solver."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

import structlog
from rdkit import Chem

from deepretro.models.hallucination_helpers import (
    HallucinationChecker,
    filter_with_checker,
    resolve_hallucination,
)
from deepretro.utils.az import run_az
from deepretro.utils.llm_helpers import Pathway
from deepretro.utils.typing import RouteNode

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

    def solve(
        self,
        smiles: str,
        visited: set[str] | None = None,
        depth: int = 0,
    ) -> tuple[dict[str, Any], bool]:
        """Recursively solve one molecule via retrosynthesis.

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

        az_solved, az_routes = run_az(smiles=smiles, az_model=self.az_model)
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
        """Call the package LLM pipeline and apply optional ML filtering.

        Parameters
        ----------
        molecule : str
            Target molecule SMILES for the LLM pipeline.

        Returns
        -------
        tuple[list[Pathway], list[str], list[float]]
            Filtered pathways, explanations, and confidence scores.
        """
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
            Solved child route nodes and a flag indicating all reactants solved.
        """
        reactants = [pathway] if isinstance(pathway, str) else list(pathway)
        children: list[dict[str, Any]] = []
        all_solved = True
        for reactant in reactants:
            child, solved = self.solve(reactant, visited, depth + 1)
            children.append(child)
            all_solved = all_solved and solved
        return children, all_solved


def canonicalize(smiles: str) -> str:
    """Return canonical SMILES, or the original string if parsing fails.

    Parameters
    ----------
    smiles : str
        Input SMILES string.

    Returns
    -------
    str
        Canonical SMILES or the original string if RDKit cannot parse it.

    Examples
    --------
    >>> canonicalize("C(O)C")
    'CCO'
    >>> canonicalize("not_a_smiles")
    'not_a_smiles'
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return smiles
    return Chem.MolToSmiles(molecule, canonical=True)


def llm_pipeline(**kwargs: Any) -> tuple[list[Pathway], list[str], list[float]]:
    """Import and call the LLM pipeline lazily.

    Keeping this import lazy lets users import ``deepretro.AutoSolver`` in
    lightweight environments that have not installed LLM provider extras yet.

    Parameters
    ----------
    **kwargs : Any
        Keyword arguments forwarded to the package LLM pipeline.

    Returns
    -------
    tuple[list[Pathway], list[str], list[float]]
        Pathways, explanations, and confidence scores from the pipeline.
    """
    from deepretro.utils.llm import llm_pipeline as package_llm_pipeline

    return package_llm_pipeline(**kwargs)


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
