"""Pipeline-compatible wrappers for package-native heuristics."""

from __future__ import annotations

from deepretro.algorithms.hallucination_checker import calculate_hallucination_score
from deepretro.utils.utils_molecule import is_valid_smiles

ALLOWED_HALLUCINATION_SEVERITIES = {"low", "medium"}


def hallucination_checker(product: str, res_smiles: list) -> tuple[int, list]:
    """Filter pathways using the package hallucination heuristic.

    Parameters
    ----------
    product : str
        SMILES string of the target product molecule.
    res_smiles : list
        Candidate pathways represented as either lists of reactant SMILES
        strings or a single reactant SMILES string.

    Returns
    -------
    tuple[int, list]
        ``(200, valid_pathways)`` when at least one valid pathway remains,
        otherwise ``(400, [])``.
    """
    valid_pathways: list[list[str]] = []

    for pathway in res_smiles:
        if isinstance(pathway, list):
            reactants_smi = ".".join(pathway)
            normalized_pathway = pathway
        else:
            reactants_smi = pathway
            normalized_pathway = [pathway]

        if not is_valid_smiles(reactants_smi):
            continue

        hallucination_report = calculate_hallucination_score(reactants_smi, product)
        if hallucination_report["severity"] in ALLOWED_HALLUCINATION_SEVERITIES:
            valid_pathways.append(normalized_pathway)

    if valid_pathways:
        return 200, valid_pathways
    return 400, []
