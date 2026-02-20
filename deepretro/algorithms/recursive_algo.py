"""Recursive function to run Prithvi on a molecule"""

from rdkit import Chem

from src.utils.az import run_az
from src.utils.job_context import logger as context_logger
from src.utils.llm import llm_pipeline


def rec_run_DeepRetro(
    molecule: str,
    job_id: str,
    llm: str = "claude-opus-4-20250514",
    az_model: str = "USPTO",
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
    visited=None,
    depth=0,
    max_depth=50,
) -> tuple[dict, bool]:
    """Recursive function to run Prithvi on a molecule

    Parameters
    ----------
    molecule : str
        Molecule SMILES
    job_id : str
        Job ID
    llm : str, optional
        LLM to be used, by default "claude-opus-4-20250514"
    az_model : str, optional
        AZ model to be used, by default "USPTO"
    stability_flag : str, optional
        Stability flag, by default "False"
    hallucination_check : str, optional
        Hallucination check, by default "False"
    use_protecting_group_feature : bool, optional
        Use protecting group feature, by default False
    visited : set, optional
        Set of molecules already processed (for cycle detection), by default None
    depth : int, optional
        Current recursion depth, by default 0
    max_depth : int, optional
        Maximum recursion depth, by default 50

    Returns
    -------
    tuple(dict, bool)
        result_dict: result of retrosynthesis.
        solved: boolean value indicating if the molecule was solved.
    """
    if visited is None:
        visited = set()
    logger = context_logger.get()

    # Canonicalize SMILES for proper cycle detection (handles different representations of same molecule)
    try:
        mol: Chem.Mol = Chem.MolFromSmiles(molecule)
        if mol:
            canonical_molecule = Chem.MolToSmiles(mol, canonical=True)
        else:
            canonical_molecule = molecule
    except Exception:
        canonical_molecule = molecule

    if depth >= max_depth:
        logger.warning(f"Max depth {max_depth} reached for {molecule}")
        return {
            "type": "mol",
            "smiles": molecule,
            "is_chemical": True,
            "in_stock": False,
            "children": [],
        }, False
    if canonical_molecule in visited:
        logger.warning(
            f"Cycle detected: {molecule} (canonical: {canonical_molecule}) already processed"
        )
        return {
            "type": "mol",
            "smiles": molecule,
            "is_chemical": True,
            "in_stock": False,
            "children": [],
        }, False
    visited.add(canonical_molecule)
    solved, result_dict = run_az(smiles=molecule, az_model=az_model)
    result_dict = result_dict[0]
    if not solved:
        logger.info(f"AZ failed for {molecule}, running LLM")
        out_pathways, out_explained, out_confidence = llm_pipeline(
            molecule=molecule,
            LLM=llm,
            stability_flag=stability_flag,
            hallucination_check=hallucination_check,
            use_protecting_group_feature=use_protecting_group_feature,
        )
        result_dict: dict[str, str | bool | list[dict]] = {
            "type": "mol",
            "smiles": molecule,
            "is_chemical": True,
            "in_stock": False,
            "children": [
                {
                    "type": "reaction",
                    "is_reaction": True,
                    "metadata": {
                        "policy_probability": out_confidence,
                    },
                    "children": [],
                }
            ],
        }
        logger.info(f"LLM returned {out_pathways}")
        logger.info(f"LLM explained {out_explained}")
        for pathway in out_pathways:
            if isinstance(pathway, list):
                temp_stat = []
                for molecule in pathway:
                    res, stat = rec_run_DeepRetro(
                        molecule=molecule,
                        job_id=job_id,
                        llm=llm,
                        az_model=az_model,
                        stability_flag=stability_flag,
                        hallucination_check=hallucination_check,
                        use_protecting_group_feature=use_protecting_group_feature,
                        visited=visited,
                        depth=depth + 1,
                        max_depth=max_depth,
                    )
                    if stat:
                        temp_stat.append(True)
                        result_dict["children"][0]["children"].append(res)
                logger.info(f"temp_stat: {temp_stat}")
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
                    depth=depth + 1,
                    max_depth=max_depth,
                )
                result_dict["children"][0]["children"].append(res)
            if solved:
                logger.info("breaking")
                break
    else:
        logger.info(f"AZ solved {molecule}")
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved


def single_run_DeepRetro(
    molecule: str,
    llm: str = "anthropic/claude-opus-4-20250514",
    az_model: str = "USPTO",
    stability_flag: str = "False",
    hallucination_check: str = "False",
    use_protecting_group_feature: bool = False,
) -> tuple[dict, bool]:
    """Single run function to run DeepRetro on a molecule

    Parameters
    ----------
    molecule : str
        Molecule SMILES
    llm : str, optional
        LLM to be used, by default "claude-opus-4-20250514"
    az_model : str, optional
        AZ model to be used, by default "USPTO"
    stability_flag : str, optional
        Stability flag, by default "False"
    hallucination_check : str, optional
        Hallucination check, by default "False"

    Returns
    -------
    tuple(dict, bool)
        result_dict: result of retrosynthesis.
        solved: boolean value indicating if the molecule was solved.
    """
    solved, result_dict = run_az(smiles=molecule, az_model=az_model)
    result_dict = result_dict[0]
    logger = context_logger.get()
    logger.info(f"AZ failed for {molecule}, running LLM")
    out_pathways, out_explained, out_confidence = llm_pipeline(
        molecule=molecule,
        LLM=llm,
        stability_flag=stability_flag,
        hallucination_check=hallucination_check,
        use_protecting_group_feature=use_protecting_group_feature,
    )
    result_dict = {
        "type": "mol",
        "smiles": molecule,
        "is_chemical": True,
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "is_reaction": True,
                "metadata": {
                    "policy_probability": out_confidence,
                },
                "children": [],
            }
        ],
    }
    logger.info(f"LLM returned {out_pathways}")
    logger.info(f"LLM explained {out_explained}")

    return result_dict, solved
