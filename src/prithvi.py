""" Module to run prithvi retrosynthesis on a molecule """
import time
import os
import rootutils
import structlog

from src.utils.parse import format_output
from src.rec_prithvi import rec_run_prithvi
from src.utils.job_context import logger as context_logger
from src.utils.custom_logging import add_job_specific_handler
from src.metadata import reagent_agent, conditions_agent, literature_agent

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

date_dir = f'{root_dir}/logs/{time.strftime("%Y-%m-%d")}'


def run_prithvi(molecule: str,
                llm="claude-opus-4-20250514",
                az_model: str = "USPTO",
                stability_flag: str = "False",
                hallucination_check: str = "False",
                use_protecting_group_feature: bool = False,
                protecting_group_mode: str = "auto",
                protecting_group_selections: dict = None,
                protection_state=None) -> dict:
    """Run prithvi services to generate retrosynthesis on a molecule.

    Parameters
    ----------
    molecule : str
        SMILE String of the molecule.
    llm : str, optional
        LLM Model, by default "claude-opus-4-20250514"
    az_model : str, optional
        AZ model to use, by default "USPTO"
    stability_flag : str, optional
        Enable stability checking, by default "False"
    hallucination_check : str, optional
        Enable hallucination checking, by default "False"
    use_protecting_group_feature : bool, optional
        Enable protecting group feature, by default False
    protecting_group_mode : str, optional
        Mode for protecting group selection: "auto" or "hitl", by default "auto"
    protecting_group_selections : dict, optional
        User-selected protecting groups keyed by site_id (HITL mode), by default None
    protection_state : ProtectionState, optional
        Tracks active protecting groups across the synthesis tree.

    Returns
    -------
    dict
        Result after running prithvi.
    """

    # Generate a unique job ID using timestamp and a random suffix
    job_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"

    job_log_file = f"{date_dir}/job_{job_id}.log"
    log = structlog.get_logger().bind(job_id=job_id)
    # Set the logger in the context variable
    token = context_logger.set(log)

    # Add job-specific handler
    handler = add_job_specific_handler(log, job_id)
    log.info(f"Starting new synthesis job {job_id} for molecule {molecule}")

    try:
        result_dict, _ = rec_run_prithvi(
            molecule=molecule,
            job_id=job_id,
            llm=llm,
            az_model=az_model,
            stability_flag=stability_flag,
            hallucination_check=hallucination_check,
            use_protecting_group_feature=use_protecting_group_feature,
            protecting_group_mode=protecting_group_mode,
            protecting_group_selections=protecting_group_selections,
            protection_state=protection_state)
        output_data = format_output(result_dict,
                                    protection_state=protection_state)
        output_data = add_metadata(output_data)
        return output_data
    finally:
        # Clean up handlers
        log._logger.removeHandler(handler)
        handler.close()
        context_logger.reset(token)


def add_metadata(output_data: dict) -> dict:
    """method to add metadata to reaction metrics

    Parameters
    ----------
    output_data : dict
        json output without metadata

    Returns
    -------
    dict
        json output with metadata
    """
    for idx, step in enumerate(output_data['steps']):
        status, reagents = reagent_agent(step['reactants'], step['products'])
        output_data['steps'][idx]['reagents'].extend(reagents)

        status, conditions = conditions_agent(step['reactants'],
                                              step['products'],
                                              step['reagents'])
        output_data['steps'][idx]['conditions'] = conditions

        status, literature = literature_agent(step['reactants'],
                                              step['products'],
                                              step['reagents'],
                                              step['conditions'])
        output_data['steps'][idx]['reactionmetrics'][0][
            'closestliterature'] = literature

    return output_data
