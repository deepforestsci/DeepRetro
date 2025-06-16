"""Implements the core recursive retrosynthesis logic for the Prithvi workflow.

This module contains the `rec_run_prithvi` function, which drives the
recursive exploration of synthetic pathways. It employs a hybrid strategy:

1.  **AiZynthFinder First:** Attempts to find a one-step retrosynthesis using
    `src.utils.az.run_az()`. If successful, the molecule is considered solved at
    this step, and its AiZynthFinder-derived precursors are used.

2.  **LLM Fallback & Recursion:** If AiZynthFinder fails to find a route for the
    current target molecule, the function queries a Large Language Model (LLM)
    via `src.utils.llm.llm_pipeline()` to propose precursor molecules (pathways).
    For each proposed precursor molecule in a pathway, `rec_run_prithvi` is called
    recursively.

3.  **Route Construction:** The function builds a tree-like dictionary (`result_dict`)
    representing the explored synthetic routes. Nodes in this tree are either
    molecules ('mol') or reactions ('reaction').

4.  **Success Propagation:** A molecule is considered 'solved' if AiZynthFinder resolves
    it directly, or if all precursors in at least one of its LLM-suggested pathways
    are recursively solved.

The process continues until molecules are resolved to known building blocks (implicitly,
when `run_az` indicates `solved=True` because the molecule is basic or in stock)
_or_ no further valid pathways are found.

The `job_id` is used to pass logging context via `src.utils.job_context`.
"""

from src.utils.llm import llm_pipeline
from src.utils.az import run_az
from src.utils.job_context import logger as context_logger


def rec_run_prithvi(molecule: str,
                    job_id: str,
                    llm: str = "claude-3-opus-20240229",
                    az_model: str = "USPTO",
                    stability_flag: str = "False",
                    hallucination_check: str = "False") -> tuple[dict, bool]:
    """Recursively attempts to find a retrosynthetic route for a molecule.

    Tries AiZynthFinder first. If it fails, uses an LLM to propose precursor
    sets (pathways) and recursively calls itself for each precursor. The first
    fully solved pathway from the LLM is accepted.

    Args:
        molecule (str):
            SMILES string of the target molecule for the current recursive step.
        job_id (str):
            Unique identifier for the overall synthesis job, used for logging context.
        llm (str, optional):
            Identifier for the LLM. Defaults to "claude-3-opus-20240229".
        az_model (str, optional):
            Identifier for the AiZynthFinder model. Defaults to "USPTO".
        stability_flag (str, optional):
            Flag ("True" or "False") for stability checks in the LLM pipeline.
            Defaults to "False".
        hallucination_check (str, optional):
            Flag ("True" or "False") for hallucination checks in the LLM pipeline.
            Defaults to "False".

    Returns:
        tuple[dict, bool]:

            -   `result_dict` (dict): A tree-like dictionary representing the explored
                synthetic pathway for the input `molecule`.
                Structure Example:
                ::

                    {
                        'type': 'mol',
                        'smiles': <molecule_smiles>,
                        'is_chemical': True,
                        'in_stock': <True if solved by AZ as basic/stock, False otherwise>,
                        'children': [ # Populated if LLM is used
                            {
                                'type': 'reaction',
                                'is_reaction': True,
                                'metadata': {'policy_probability': <LLM confidence>},
                                'children': [ # Precursor molecule nodes (recursive results)
                                    { 'type': 'mol', 'smiles': <precursor1_smiles>, ... },
                                    { 'type': 'mol', 'smiles': <precursor2_smiles>, ... }
                                ]
                            }
                        ]
                    }

                If AiZynthFinder solves the molecule, `result_dict` is the direct,
                formatted output from `run_az()` (which might have a different but
                compatible structure, often simpler as it's a one-step result).
            -   `solved` (bool): True if the `molecule` was successfully retrosynthesized
                either by AiZynthFinder directly, or by finding a complete recursive
                pathway through LLM-suggested precursors. False otherwise.
    """
    solved, result_dict = run_az(smiles=molecule, az_model=az_model)
    result_dict = result_dict[0]
    logger = context_logger.get()
    if not solved:
        logger.info(f"AZ failed for {molecule}, running LLM")
        out_pathways, out_explained, out_confidence = llm_pipeline(
            molecule=molecule,
            LLM=llm,
            stability_flag=stability_flag,
            hallucination_check=hallucination_check)
        result_dict = {
            'type':
            'mol',
            'smiles':
            molecule,
            # 'confidence': out_confidence,
            "is_chemical":
            True,
            "in_stock":
            False,
            'children': [{
                "type": "reaction",
                "is_reaction": True,
                "metadata": {
                    "policy_probability": out_confidence,
                },
                "children": []
            }]
        }
        logger.info(f"LLM returned {out_pathways}")
        logger.info(f"LLM explained {out_explained}")
        for pathway in out_pathways:
            if isinstance(pathway, list):
                temp_stat = []
                for mol in pathway:
                    res, stat = rec_run_prithvi(
                        molecule=mol,
                        job_id=job_id,
                        llm=llm,
                        az_model=az_model,
                        stability_flag=stability_flag,
                        hallucination_check=hallucination_check)
                    if stat:
                        temp_stat.append(True)
                        result_dict['children'][0]['children'].append(res)
                logger.info(f"temp_stat: {temp_stat}")
                if all(temp_stat):
                    solved = True
            else:
                res, solved = rec_run_prithvi(
                    molecule=pathway,
                    job_id=job_id,
                    llm=llm,
                    az_model=az_model,
                    stability_flag=stability_flag,
                    hallucination_check=hallucination_check)
                result_dict['children'][0]['children'].append(res)
            if solved:
                logger.info('breaking')
                break
    else:
        logger.info(f"AZ solved {molecule}")
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved
