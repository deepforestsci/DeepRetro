""" Recursive function to run Prithvi on a molecule """
import os
from src.utils.llm import llm_pipeline
from src.utils.az import run_az
from src.utils.job_context import logger as context_logger

ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING",
                                    "true").lower() == "false" else True

def log_message(message: str, logger=None):
    """Log the message

    Parameters
    ----------
    message : str
        The message to be logged
    logger : _type_, optional
        The logger object, by default None

    Returns
    -------
    None
    """
    if logger is not None:
        logger.info(message)
    else:
        print(message)

def rec_run_prithvi(molecule: str,
                    job_id: str,
                    llm: str = "claude-3-opus-20240229",
                    az_model: str = "USPTO",
                    stability_flag: str = "False",
                    hallucination_check: str = "False") -> tuple[dict, bool]:
    """Recursive function to run Prithvi on a molecule

    Parameters
    ----------
    molecule : str
        Molecule SMILES
    job_id : str
        Job ID
    llm : str, optional
        LLM to be used, by default "claude-3-opus-20240229"

    Returns
    -------
    tuple(dict, bool)
        result_dict: result of retrosynthesis.
        solved: boolean value indicating if the molecule was solved.
    """
    solved, result_dict = run_az(smiles=molecule, az_model=az_model)
    result_dict = result_dict[0]
    logger = context_logger.get() if ENABLE_LOGGING else None
    if not solved:
        log_message(f"AZ failed for {molecule}, running LLM", logger)
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
        log_message(f"LLM returned {out_pathways}", logger)
        log_message(f"LLM explained {out_explained}", logger)
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
                log_message(f"temp_stat: {temp_stat}", logger)
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
                log_message('breaking', logger)
                break
    else:
        log_message(f"AZ solved {molecule}", logger)
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved
