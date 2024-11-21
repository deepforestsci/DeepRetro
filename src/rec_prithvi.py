""" Recursive function to run Prithvi on a molecule """

from src.utils.llm import llm_pipeline
from src.utils.az import run_az
from src.utils.job_context import logger as context_logger


def rec_run_prithvi(molecule, job_id):
    solved, result_dict = run_az(molecule)
    result_dict = result_dict[0]
    logger = context_logger.get()
    if not solved:
        logger.info(f"AZ failed for {molecule}, running LLM")
        out_pathways, out_explained, out_confidence = llm_pipeline(molecule)
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
                    res, stat = rec_run_prithvi(mol, job_id)
                    if stat:
                        temp_stat.append(True)
                        result_dict['children'][0]['children'].append(res)
                logger.info(f"temp_stat: {temp_stat}")
                if all(temp_stat):
                    solved = True
            else:
                res, solved = rec_run_prithvi(pathway, job_id)
                result_dict['children'][0]['children'].append(res)
            if solved:
                logger.info('breaking')
                break
    else:
        logger.info(f"AZ solved {molecule}")
    # print(f"Solved : {solved}, Returning {result_dict}")
    return result_dict, solved
