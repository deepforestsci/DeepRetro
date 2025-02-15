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


def run_prithvi(molecule:str,
                llm="claude-3-opus-20240229"):
    """Run prithvi services to generate retrosynthesis on a molecule.

    Parameters
    ----------
    molecule : str
        SMILE String of the molecule. 
    llm : str, optional
        LLM Model, by default "claude-3-opus-20240229"

    Returns
    ------- #TODO: Complete docstring properly.
    _type_
        _description_
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
        result_dict, solved = rec_run_prithvi(molecule, job_id, llm)
        output_data = format_output(result_dict)
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
