"""Single-molecule retrosynthesis runner.

1. **AiZynthFinder** attempts template-based retrosynthesis.
2. If AZ fails, the **LLM** proposes retrosynthetic steps.
3. Pathways are validated (SMILES validity, stability, hallucination).
4. Recurse on each reactant until all leaves are purchasable.
5. Flatten the tree into a step-by-step synthesis plan.

Examples
--------
>>> from deepretro.algorithms.autosolve import autosolve  # doctest: +SKIP
>>> result = autosolve("CC(=O)Oc1ccccc1C(=O)O")           # doctest: +SKIP
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def autosolve(
    smiles: str,
    *,
    llm: str = "anthropic/claude-opus-4-6:adv:think",
    az_model: str = "Pistachio_100+",
    stability_check: bool = True,
    hallucination_mode: str = "heuristic",
    hallucination_classifier: Any = None,
    use_protecting_group_feature: bool = False,
    return_image: bool = False,
) -> dict[str, Any]:
    """Run retrosynthesis on a single molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.
    llm : str
        LLM model identifier.  ``"provider/model:adv:think"``
    az_model : str
        AiZynthFinder model variant (``"USPTO"``, ``"Pistachio_100+"``).
    stability_check : bool
        Discard unstable pathways (ring strain, carbocations, etc.).
    hallucination_mode : str
        ``"heuristic"`` | ``"ml"`` | ``"none"``.
    hallucination_classifier : HallucinationClassifier or str or Path or None
        Required when *hallucination_mode* is ``"ml"``.
    use_protecting_group_feature : bool
        Include protecting-group context in the LLM prompt.
    return_image : bool
        If ``True``, add an ``"image"`` key to the result containing a
        ``PIL.Image.Image`` of the synthesis pathway.

    Returns
    -------
    dict[str, Any]
        ``{"steps": [...], "dependencies": {...}}``.
        When *return_image* is ``True``, also contains ``"image"``.
    """
    import os
    import time

    import structlog
    from deepretro.logging import logger as context_logger
    from deepretro.utils.parse import format_output

    hallucination_checker = resolve_hallucination(
        hallucination_mode, hallucination_classifier,
    )

    solver = AutoSolver(
        llm=llm,
        az_model=az_model,
        stability_flag=str(stability_check),
        hallucination_checker=hallucination_checker,
        use_protecting_group_feature=use_protecting_group_feature,
    )

    job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
    log = structlog.get_logger().bind(job_id=job_id)
    token = context_logger.set(log)
    try:
        result_dict, _ = solver.solve(smiles)
        output = format_output(result_dict)
    finally:
        context_logger.reset(token)

    if return_image:
        from deepretro.utils.visualize import visualize_pathway
        output["image"] = visualize_pathway(output)

    return output


class AutoSolver:
    """Holds shared config"""

    def __init__(
        self,
        llm: str,
        az_model: str,
        stability_flag: str,
        hallucination_checker,
        use_protecting_group_feature: bool,
        max_depth: int = 50,
    ):
        self.llm = llm
        self.az_model = az_model
        self.stability_flag = stability_flag
        self.hallucination_checker = hallucination_checker
        self.use_protecting_group_feature = use_protecting_group_feature
        self.max_depth = max_depth

    def solve(
        self,
        molecule: str,
        visited: set | None = None,
        depth: int = 0,
    ) -> tuple[dict, bool]:
        """Recursively solve a molecule via AZ, falling back to LLM.

        Returns a nested mol/reaction tree and a *solved* flag.
        """
        from rdkit import Chem
        from deepretro.algorithms.llm import llm_pipeline
        from deepretro.utils.az import run_az
        from deepretro.logging import logger as context_logger

        logger = context_logger.get()

        if visited is None:
            visited = set()

        try:
            mol = Chem.MolFromSmiles(molecule)
            canonical = Chem.MolToSmiles(mol, canonical=True) if mol else molecule
        except Exception:
            canonical = molecule

        if depth >= self.max_depth:
            logger.warning(f"Max depth {self.max_depth} reached for {molecule}")
            return unsolved_leaf(molecule), False
        if canonical in visited:
            logger.warning(f"Cycle detected: {molecule} (canonical: {canonical})")
            return unsolved_leaf(molecule), False

        visited.add(canonical)

        solved, az_results = run_az(smiles=molecule, az_model=self.az_model)
        result_dict = az_results[0]
        if solved:
            logger.info(f"AZ solved {molecule}")
            return result_dict, True

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

        for step in pathways:
            reactants = step if isinstance(step, list) else [step]
            all_solved = True
            for smi in reactants:
                res, stat = self.solve(smi, visited, depth + 1)
                if stat:
                    children.append(res)
                else:
                    all_solved = False
            solved = all_solved
            if solved:
                break

        return result_dict, solved


def unsolved_leaf(smiles: str) -> dict:
    """Terminal node for molecules that hit max depth or cycle detection."""
    return {
        'type': 'mol', 'smiles': smiles,
        'is_chemical': True, 'in_stock': False, 'children': [],
    }


def resolve_hallucination(mode: str, classifier: Any):
    """Convert user-facing mode to a single checker callable (or None).

    Returns None (skip checking), or a callable with signature
    ``(product: str, pathways: list) -> (int, list)``.
    """
    if mode == "none":
        return None
    if mode == "heuristic":
        from deepretro.algorithms.pipeline_checks import hallucination_checker
        return hallucination_checker
    if mode == "ml":
        from deepretro.models.hallucination_classifier import HallucinationClassifier
        from deepretro.models.hallucination_helpers import build_ml_checker

        if isinstance(classifier, (str, Path)):
            clf = HallucinationClassifier()
            clf.load(str(classifier))
        elif hasattr(classifier, "predict_single"):
            clf = classifier
        else:
            raise ValueError(
                f"hallucination_mode='ml' requires a HallucinationClassifier "
                f"or path to saved model — got {type(classifier)}"
            )
        return build_ml_checker(clf)
    raise ValueError(
        f"hallucination_mode must be 'heuristic', 'ml', or 'none' — got {mode!r}"
    )
