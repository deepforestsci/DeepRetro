"""Recursive retrosynthesis using DeepRetro.

Attempts AiZynthFinder (template-based) first; falls back to an LLM pipeline
when the molecule cannot be solved. Recurses on proposed reactants until
all leaves are commercially available or depth/cycle limits are hit.
"""

from __future__ import annotations

from typing import TypedDict, cast

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from deepretro.utils.az import run_az
from src.utils.job_context import logger as context_logger
from deepretro.utils.llm import llm_pipeline


def _get_logger():
    """Return the context logger if set, else None. Safe for doctest/standalone use."""
    return context_logger.get(None)


def _log(level: str, message: str, logger=None):
    """Log message at level if logger is set, else no-op."""
    if logger is not None:
        getattr(logger, level)(message)


class ReactionMetadata(TypedDict):
    """Metadata attached to a reaction node.

    Attributes
    ----------
    policy_probability : list[float]
        Confidence scores for each proposed precursor in the reaction.
    """

    policy_probability: list[float]


class ReactionNode(TypedDict):
    """A single retrosynthetic reaction step.

    Attributes
    ----------
    type : str
        Node type, typically ``"reaction"``.
    is_reaction : bool
        Always ``True`` for reaction nodes.
    metadata : ReactionMetadata
        Confidence and policy information.
    children : list[MolNode]
        Precursor molecule nodes produced by this reaction.
    """

    type: str
    is_reaction: bool
    metadata: ReactionMetadata
    children: list[MolNode]


class MolNode(TypedDict):
    """A molecule in the retrosynthesis tree.

    Attributes
    ----------
    type : str
        Node type, typically ``"mol"``.
    smiles : str
        SMILES string of the molecule.
    is_chemical : bool
        Whether the molecule is a valid chemical structure.
    in_stock : bool
        Whether the molecule is commercially available.
    children : list[ReactionNode]
        Retrosynthetic reactions that produce this molecule.
    """

    type: str
    smiles: str
    is_chemical: bool
    in_stock: bool
    children: list[ReactionNode]


ResultPair = tuple[MolNode, bool]
"""Type alias: ``(result_tree, solved)`` - retrosynthesis sub-tree and solved status."""


def _make_mol_node(
    smiles: str,
    in_stock: bool = False,
    children: list[ReactionNode] | None = None,
) -> MolNode:
    """Build a molecule node for the retrosynthesis tree.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    in_stock : bool, optional
        Whether the molecule is in stock, by default False.
    children : list[ReactionNode] | None, optional
        Children of the molecule node, by default None.

    Returns
    -------
    MolNode
        Molecule node dict with type, smiles, is_chemical, in_stock, children.

    Examples
    --------
    >>> from deepretro.algorithms.recursive_algo import _make_mol_node
    >>> node = _make_mol_node("CCO")
    >>> node["smiles"]
    'CCO'
    >>> node["type"]
    'mol'
    >>> node["in_stock"]
    False
    >>> node = _make_mol_node("CC", in_stock=True)
    >>> node["in_stock"]
    True
    """
    return {
        "type": "mol",
        "smiles": smiles,
        "is_chemical": True,
        "in_stock": in_stock,
        "children": children if children is not None else [],
    }


def _make_reaction_node(
    confidence: list[float],
    children: list[MolNode] | None = None,
) -> ReactionNode:
    """Build a reaction node for the retrosynthesis tree.

    Parameters
    ----------
    confidence : list[float]
        Confidence scores for the reaction.
    children : list[MolNode] | None, optional
        Children of the reaction node, by default None.

    Returns
    -------
    ReactionNode
        Reaction node dict with type, is_reaction, metadata, children.

    Examples
    --------
    >>> from deepretro.algorithms.recursive_algo import _make_reaction_node
    >>> node = _make_reaction_node([0.9, 0.8])
    >>> node["type"]
    'reaction'
    >>> node["metadata"]["policy_probability"]
    [0.9, 0.8]
    >>> node["is_reaction"]
    True
    """
    return {
        "type": "reaction",
        "is_reaction": True,
        "metadata": {"policy_probability": confidence},
        "children": children if children is not None else [],
    }


def _canonicalize(smiles: str) -> str:
    """Return canonical SMILES, falling back to the original on failure.

    Parameters
    ----------
    smiles : str
        SMILES string to canonicalize.

    Returns
    -------
    str
        Canonical SMILES string, or original on parse failure.

    Examples
    --------
    >>> from deepretro.algorithms.recursive_algo import _canonicalize
    >>> _canonicalize("c1ccccc1")
    'c1ccccc1'
    >>> _canonicalize("CCO")
    'CCO'
    >>> _canonicalize("invalid")
    'invalid'
    """
    try:
        mol: Mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        pass
    return smiles


def rec_run_DeepRetro(
    molecule: str,
    job_id: str,
    llm: str = "claude-opus-4-6",
    az_model: str = "USPTO",
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
    visited: set[str] | None = None,
    depth: int = 0,
    max_depth: int = 50,
) -> ResultPair:
    """Recursively decompose a molecule via retrosynthesis.

    Attempts AiZynthFinder first; falls back to an LLM pipeline when AZ
    cannot solve the molecule.  Recurses on each proposed reactant until
    all leaves are commercially available or limits are hit.

    Parameters
    ----------
    molecule : str
        Target molecule as a SMILES string.
    job_id : str
        Unique identifier for the current job.
    llm : str, optional
        LLM model name, by default ``"claude-opus-4-6"``.
    az_model : str, optional
        AiZynthFinder model variant, by default ``"USPTO"``.
    stability_flag : str, optional
        Enable stability checking, by default ``"False"``.
    hallucination_check : str, optional
        Enable hallucination checking, by default ``"False"``.
    use_protecting_group_feature : bool, optional
        Enable protecting-group suggestions, by default ``False``.
    visited : set[str] | None, optional
        Canonical SMILES already processed (cycle detection).
    depth : int, optional
        Current recursion depth, by default ``0``.
    max_depth : int, optional
        Maximum recursion depth, by default ``50``.

    Returns
    -------
    ResultPair
        ``(result_dict, solved)`` — the retrosynthesis sub-tree rooted at
        *molecule* and whether a complete route was found.

    Examples
    --------
    >>> from deepretro.algorithms.recursive_algo import rec_run_DeepRetro
    >>> from src.cache import clear_cache_for_molecule
    >>> clear_cache_for_molecule("CCO")
    >>> tree, solved = rec_run_DeepRetro("CCO", "job-1", llm='claude-haiku-4-5-20251001') # doctest: +SKIP
    >>> isinstance(tree, dict) and "smiles" in tree # doctest: +SKIP
    True
    >>> tree["smiles"] # doctest: +SKIP
    'CCO'
    """
    if visited is None:
        visited = set()
    logger = _get_logger()

    canonical = _canonicalize(molecule)

    if depth >= max_depth:
        _log("warning", f"Max depth {max_depth} reached for {molecule}", logger)
        return _make_mol_node(molecule), False

    if canonical in visited:
        _log(
            "warning",
            f"Cycle detected: {molecule} (canonical: {canonical}) already processed",
            logger,
        )
        return _make_mol_node(molecule), False

    visited.add(canonical)

    solved, az_results = run_az(smiles=molecule, az_model=az_model)
    result_dict = cast(MolNode, az_results[0])

    if solved:
        _log("info", f"AZ solved {molecule}", logger)
        return result_dict, bool(solved)

    _log("info", f"AZ failed for {molecule}, running LLM", logger)
    out_pathways, out_explained, out_confidence = llm_pipeline(
        molecule=molecule,
        model=llm,
        stability_check=stability_flag,
        hallucination_check=hallucination_check,
        use_protecting_group_feature=use_protecting_group_feature,
    )

    reaction = _make_reaction_node(out_confidence)
    result_dict = _make_mol_node(molecule, children=[reaction])

    _log("info", f"LLM returned {out_pathways}", logger)
    _log("info", f"LLM explained {out_explained}", logger)

    next_depth = depth + 1

    for pathway in out_pathways:
        if isinstance(pathway, list):
            temp_stat: list[bool] = []
            for reactant in pathway:
                res, stat = rec_run_DeepRetro(
                    molecule=reactant,
                    job_id=job_id,
                    llm=llm,
                    az_model=az_model,
                    stability_flag=stability_flag,
                    hallucination_check=hallucination_check,
                    use_protecting_group_feature=use_protecting_group_feature,
                    visited=visited,
                    depth=next_depth,
                    max_depth=max_depth,
                )
                if stat:
                    temp_stat.append(True)
                    reaction["children"].append(res)
            _log("info", f"temp_stat: {temp_stat}", logger)
            if all(temp_stat):
                solved = True
        else:
            res, solved = rec_run_DeepRetro(
                molecule=pathway,
                job_id=job_id,
                llm=llm,
                az_model=az_model,
                stability_flag=stability_flag,
                hallucination_check=hallucination_check,
                use_protecting_group_feature=use_protecting_group_feature,
                visited=visited,
                depth=next_depth,
                max_depth=max_depth,
            )
            reaction["children"].append(res)

        if solved:
            _log("info", "breaking", logger)
            break

    return result_dict, bool(solved)


def single_run_DeepRetro(
    molecule: str,
    llm: str = "anthropic/claude-opus-4-6",
    az_model: str = "USPTO",
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
) -> ResultPair:
    """Run a single (non-recursive) retrosynthesis step on a molecule.

    Always runs both AiZynthFinder and the LLM pipeline.  Returns the
    AZ solved status alongside the LLM-produced result tree.

    Parameters
    ----------
    molecule : str
        Target molecule as a SMILES string.
    llm : str, optional
        LLM model name, by default ``"anthropic/claude-opus-4-6"``.
    az_model : str, optional
        AiZynthFinder model variant, by default ``"USPTO"``.
    stability_flag : str, optional
        Enable stability checking, by default ``"False"``.
    hallucination_check : str, optional
        Enable hallucination checking, by default ``"False"``.
    use_protecting_group_feature : bool, optional
        Enable protecting-group suggestions, by default ``False``.

    Returns
    -------
    ResultPair
        ``(result_dict, solved)`` — the LLM-produced retrosynthesis node
        and the AZ solved status.

    Examples
    --------
    >>> from deepretro.algorithms.recursive_algo import single_run_DeepRetro
    >>> from src.cache import clear_cache_for_molecule
    >>> clear_cache_for_molecule("CCO")
    >>> tree, solved = single_run_DeepRetro("CCO", llm="claude-haiku-4-5-20251001") # doctest: +SKIP
    >>> isinstance(tree, dict) and "children" in tree # doctest: +SKIP
    True
    >>> tree["smiles"] # doctest: +SKIP
    'CCO'
    """
    logger = _get_logger()

    solved, _az_results = run_az(smiles=molecule, az_model=az_model)

    _log("info", f"AZ failed for {molecule}, running LLM", logger)
    out_pathways, out_explained, out_confidence = llm_pipeline(
        molecule=molecule,
        model=llm,
        stability_check=stability_flag,
        hallucination_check=hallucination_check,
        use_protecting_group_feature=use_protecting_group_feature,
    )

    reaction = _make_reaction_node(out_confidence)
    result_dict = _make_mol_node(molecule, children=[reaction])

    _log("info", f"LLM returned {out_pathways}", logger)
    _log("info", f"LLM explained {out_explained}", logger)

    return result_dict, bool(solved)
