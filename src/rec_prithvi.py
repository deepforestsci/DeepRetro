""" Recursive function to run Prithvi on a molecule """

from src.utils.llm import llm_pipeline, agentic_llm_pipeline
from src.utils.az import run_az
from src.utils.job_context import logger as context_logger


def rec_run_prithvi(
        molecule: str,
        job_id: str,
        llm: str = "claude-opus-4-20250514",
        az_model: str = "USPTO",
        stability_flag: str = "False",
        hallucination_check: str = "False",
        use_protecting_group_feature: bool = False,
        use_agentic_pipeline: bool = True,
        max_tool_iterations: int = 5) -> tuple[dict, bool]:
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
        Whether to use protecting group feature, by default False
    use_agentic_pipeline : bool, optional
        Whether to use agentic pipeline with tools, by default True
    max_tool_iterations : int, optional
        Maximum number of tool iterations for agentic pipeline, by default 5

    Returns
    -------
    tuple(dict, bool)
        result_dict: result of retrosynthesis.
        solved: boolean value indicating if the molecule was solved.
    """
    solved, result_dict = run_az(smiles=molecule, az_model=az_model)
    result_dict = result_dict[0]
    logger = context_logger.get()
    if not solved:
        logger.info(f"AZ failed for {molecule}, running {'Agentic' if use_agentic_pipeline else 'Standard'} LLM")
        
        # Use agentic or standard pipeline based on flag
        if use_agentic_pipeline:
            out_pathways, out_explained, out_confidence, tool_stats = agentic_llm_pipeline(
                molecule=molecule,
                LLM=llm,
                stability_flag=stability_flag,
                hallucination_check=hallucination_check,
                use_protecting_group_feature=use_protecting_group_feature,
                max_tool_iterations=max_tool_iterations,
                trace_id=job_id
            )
            logger.info(f"Tool usage statistics: {tool_stats}")
        else:
            out_pathways, out_explained, out_confidence = llm_pipeline(
                molecule=molecule,
                LLM=llm,
                stability_flag=stability_flag,
                hallucination_check=hallucination_check,
                use_protecting_group_feature=use_protecting_group_feature)
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
                        hallucination_check=hallucination_check,
                        use_protecting_group_feature=use_protecting_group_feature,
                        use_agentic_pipeline=use_agentic_pipeline,
                        max_tool_iterations=max_tool_iterations)
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
                    hallucination_check=hallucination_check,
                    use_protecting_group_feature=use_protecting_group_feature,
                    use_agentic_pipeline=use_agentic_pipeline,
                    max_tool_iterations=max_tool_iterations)
                result_dict['children'][0]['children'].append(res)
            if solved:
                logger.info('breaking')
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
        use_agentic_pipeline: bool = True,
        max_tool_iterations: int = 5) -> tuple[dict, bool]:
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
    use_protecting_group_feature : bool, optional
        Whether to use protecting group feature, by default False
    use_agentic_pipeline : bool, optional
        Whether to use agentic pipeline with tools, by default True
    max_tool_iterations : int, optional
        Maximum number of tool iterations for agentic pipeline, by default 5

    Returns
    -------
    tuple(dict, bool)
        result_dict: result of retrosynthesis.
        solved: boolean value indicating if the molecule was solved.
    """
    solved, result_dict = run_az(smiles=molecule, az_model=az_model)
    result_dict = result_dict[0]
    logger = context_logger.get()
    logger.info(f"AZ failed for {molecule}, running {'Agentic' if use_agentic_pipeline else 'Standard'} LLM")
    
    # Use agentic or standard pipeline based on flag
    if use_agentic_pipeline:
        out_pathways, out_explained, out_confidence, tool_stats = agentic_llm_pipeline(
            molecule=molecule,
            LLM=llm,
            stability_flag=stability_flag,
            hallucination_check=hallucination_check,
            use_protecting_group_feature=use_protecting_group_feature,
            max_tool_iterations=max_tool_iterations
        )
        logger.info(f"Tool usage statistics: {tool_stats}")
    else:
        out_pathways, out_explained, out_confidence = llm_pipeline(
            molecule=molecule,
            LLM=llm,
            stability_flag=stability_flag,
            hallucination_check=hallucination_check,
            use_protecting_group_feature=use_protecting_group_feature)
    result_dict = {
        'type':
        'mol',
        'smiles':
        molecule,
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

    return result_dict, solved
