"""Single-molecule retrosynthesis solver.

1. **AiZynthFinder** attempts template-based retrosynthesis.
2. If AZ fails, the **LLM** proposes retrosynthetic steps.
3. Pathways are validated (SMILES validity, stability, hallucination).
4. Recurse on each reactant until all leaves are purchasable.
5. Flatten the tree into a step-by-step synthesis plan.

Examples
--------
>>> from deepretro.algorithms.autosolve import AutoSolver  # doctest: +SKIP
>>> solver = AutoSolver()                                   # doctest: +SKIP
>>> result = solver.solve("CC(=O)Oc1ccccc1C(=O)O")         # doctest: +SKIP
"""

from __future__ import annotations

import os
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import structlog
from rdkit import Chem

from deepretro.algorithms.llm import llm_pipeline
from deepretro.logging import logger as context_logger
from deepretro.models.hallucination_helpers import resolve_hallucination
from deepretro.utils.az import run_az
from deepretro.utils.parse import format_output


class AutoSolver:
    """Retrosynthesis solver that wraps AiZynthFinder + LLM fallback.

    Create once, call :meth:`solve` for each molecule.  Configuration
    is stored on the instance so repeated calls avoid redundant setup
    (e.g. reloading an ML hallucination model).

    Parameters
    ----------
    llm : str
        LLM model identifier.  ``"provider/model:adv:think"``
    az_model : str
        AiZynthFinder model variant (``"USPTO"``, ``"Pistachio_100+"``).
    stability_check : bool
        Discard unstable pathways (ring strain, carbocations, etc.).
    hallucination_mode : str
        ``"heuristic"`` | ``"ml"`` | ``"none"``.
    hallucination_classifier : str or Path or None
        Path to a saved ML model directory.
        Required when *hallucination_mode* is ``"ml"``.
    use_protecting_group_feature : bool
        Include protecting-group context in the LLM prompt.
    max_depth : int
        Maximum recursion depth for the solver.

    Examples
    --------
    >>> solver = AutoSolver()                           # doctest: +SKIP
    >>> result = solver.solve("CCO")                    # doctest: +SKIP
    >>> solver = AutoSolver(hallucination_mode="ml",    # doctest: +SKIP
    ...     hallucination_classifier="model_out/")
    >>> for smi in smiles_list:                         # doctest: +SKIP
    ...     result = solver.solve(smi)
    """

    def __init__(
        self,
        llm: str = "anthropic/claude-opus-4-6:adv:think",
        az_model: str = "Pistachio_100+",
        stability_check: bool = True,
        hallucination_mode: str = "heuristic",
        hallucination_classifier: str | Path | None = None,
        use_protecting_group_feature: bool = False,
        max_depth: int = 50,
    ):
        self.llm = llm
        self.az_model = az_model
        self.stability_flag = str(stability_check)
        self.hallucination_checker: Callable | None = resolve_hallucination(
            hallucination_mode, hallucination_classifier,
        )
        self.use_protecting_group_feature = use_protecting_group_feature
        self.max_depth = max_depth

    def solve(
        self,
        smiles: str,
        *,
        return_image: bool = False,
    ) -> dict[str, Any]:
        """Run retrosynthesis on a single molecule.

        Parameters
        ----------
        smiles : str
            SMILES string of the target molecule.
        return_image : bool
            If ``True``, add an ``"image"`` key containing a
            ``PIL.Image.Image`` of the synthesis pathway.

        Returns
        -------
        dict[str, Any]
            ``{"steps": [...], "dependencies": {...}}``.
            When *return_image* is ``True``, also contains ``"image"``.

        Examples
        --------
        >>> solver = AutoSolver(hallucination_mode="none")  # doctest: +SKIP
        >>> result = solver.solve("CCO")                    # doctest: +SKIP
        >>> result.keys()                                   # doctest: +SKIP
        dict_keys(['steps', 'dependencies'])
        >>> result = solver.solve("CCO", return_image=True) # doctest: +SKIP
        >>> "image" in result                               # doctest: +SKIP
        True
        """
        job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
        log = structlog.get_logger().bind(job_id=job_id)
        token = context_logger.set(log)
        try:
            result_dict, _ = self.recurse(smiles)
            output = format_output(result_dict)
        finally:
            context_logger.reset(token)

        if return_image:
            from deepretro.utils.visualize import visualize_pathway
            output["image"] = visualize_pathway(output)

        return output

    def recurse(
        self,
        molecule: str,
        visited: set | None = None,
        depth: int = 0,
    ) -> tuple[dict, bool]:
        """Recursively solve a molecule via AZ, falling back to LLM.

        For each molecule the solver first tries AiZynthFinder (template-based).
        If AZ finds no route, the LLM proposes candidate reaction steps.  Each
        reactant in the LLM's response is then solved recursively until every
        leaf is purchasable or the depth limit is reached.

        Parameters
        ----------
        molecule : str
            SMILES string of the molecule to solve.
        visited : set or None
            Canonical SMILES already seen on this branch (cycle detection).
            Initialised automatically on the first call.
        depth : int
            Current recursion depth.  Stops at ``self.max_depth``.

        Returns
        -------
        tuple[dict, bool]
            Nested mol/reaction tree and a *solved* flag.

        Examples
        --------
        >>> solver = AutoSolver(hallucination_mode="none")  # doctest: +SKIP
        >>> tree, solved = solver.recurse("CCO")            # doctest: +SKIP
        """
        logger = context_logger.get()
        if visited is None:
            visited = set()

        canonical = canonicalize(molecule)

        if depth >= self.max_depth:
            logger.warning(f"Max depth {self.max_depth} reached for {molecule}")
            return unsolved_leaf(molecule), False
        if canonical in visited:
            logger.warning(f"Cycle detected: {molecule} (canonical: {canonical})")
            return unsolved_leaf(molecule), False

        descendant_visited = set(visited)
        descendant_visited.add(canonical)

        solved, az_results = run_az(smiles=molecule, az_model=self.az_model)
        if solved:
            logger.info(f"AZ solved {molecule}")
            return az_results[0], True

        logger.info(f"AZ failed for {molecule}, running LLM")
        pathways, explained, confidence = llm_pipeline(
            molecule=molecule,
            LLM=self.llm,
            stability_flag=self.stability_flag,
            hallucination_checker=self.hallucination_checker,
            use_protecting_group_feature=self.use_protecting_group_feature,
        )
        logger.info(f"LLM returned {pathways}")
        logger.info(f"LLM explained {explained}")

        result_dict = {
            'type': 'mol',
            'smiles': molecule,
            'is_chemical': True,
            'in_stock': False,
            'children': [{
                'type': 'reaction',
                'is_reaction': True,
                'metadata': {'policy_probability': confidence},
                'children': [],
            }],
        }
        children = result_dict['children'][0]['children']
        all_solved = False

        for step in pathways:
            reactants = step if isinstance(step, list) else [step]
            candidate_children = []
            candidate_solved = True
            for smi in reactants:
                res, stat = self.recurse(smi, descendant_visited.copy(), depth + 1)
                if stat:
                    candidate_children.append(res)
                else:
                    candidate_solved = False
            if candidate_solved:
                children.extend(candidate_children)
                all_solved = True
                break

        return result_dict, all_solved


def canonicalize(smiles: str) -> str:
    """Return canonical SMILES, or the original string on failure.

    Parameters
    ----------
    smiles : str
        Input SMILES string (may be non-canonical or invalid).

    Returns
    -------
    str
        Canonical SMILES if RDKit can parse the input, otherwise the
        original string unchanged.

    Examples
    --------
    >>> canonicalize("C(O)C")
    'CCO'
    >>> canonicalize("INVALID")
    'INVALID'
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True) if mol else smiles
    except Exception:
        return smiles



def unsolved_leaf(smiles: str) -> dict:
    """Build a terminal tree node for unsolved molecules.

    Used when a molecule hits the max recursion depth or is part of
    a cycle, so the solver cannot continue.

    Parameters
    ----------
    smiles : str
        SMILES string for the unsolved molecule.

    Returns
    -------
    dict
        Leaf node with ``"type": "mol"`` and no children.

    Examples
    --------
    >>> leaf = unsolved_leaf("CCO")
    >>> leaf["children"]
    []
    >>> leaf["in_stock"]
    False
    """
    return {
        'type': 'mol', 'smiles': smiles,
        'is_chemical': True, 'in_stock': False, 'children': [],
    }
