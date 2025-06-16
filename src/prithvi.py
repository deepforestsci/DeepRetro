""" Orchestrates the Prithvi retrosynthesis workflow for a target molecule.

This module provides the main `run_prithvi` function, which coordinates the
end-to-end retrosynthesis process using the "Prithvi" methodology. This involves:

1.  **Job-Specific Logging:** Setting up a unique logger for each synthesis job,
    with outputs directed to a dated log file (`logs/<date>/job_<job_id>.log`).
    A context logger is used to pass this job-specific logger through the call stack.
2.  **Recursive Retrosynthesis:** Calling `src.rec_prithvi.rec_run_prithvi()` to
    perform the core, potentially recursive, retrosynthetic analysis.
3.  **Output Formatting:** Structuring the raw results from the recursive step into a
    standardized format using `src.utils.parse.format_output()`.
4.  **Metadata Enrichment:** Augmenting each reaction step in the pathway with
    predicted reagents, reaction conditions, and literature references using
    agents from `src.metadata` via the `add_metadata()` helper function.

The final output is a dictionary containing the complete retrosynthetic pathway
with all associated predictions and metadata.
"""
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
                llm="claude-3-opus-20240229",
                az_model: str = "USPTO",
                stability_flag: str = "False",
                hallucination_check: str = "False") -> dict:
    """Runs the Prithvi retrosynthesis workflow for a given molecule.

    This function orchestrates the retrosynthesis process, including setting up
    job-specific logging, calling the recursive Prithvi algorithm, formatting
    the output, and enriching it with predicted metadata.

    Args:
        molecule (str):
            The SMILES string of the target molecule for retrosynthesis.
        llm (str, optional):
            Identifier for the Large Language Model to be used. Defaults to
            "claude-3-opus-20240229". Passed to `rec_run_prithvi`.
        az_model (str, optional):
            Identifier for the AiZynthFinder model/configuration. Defaults to
            "USPTO". Passed to `rec_run_prithvi`.
        stability_flag (str, optional):
            Flag ("True" or "False") to enable/disable stability checks.
            Defaults to "False". Passed to `rec_run_prithvi`.
        hallucination_check (str, optional):
            Flag ("True" or "False") to enable/disable hallucination checks.
            Defaults to "False". Passed to `rec_run_prithvi`.

    Returns:
        dict:
            A dictionary representing the final retrosynthesis pathway and metadata.
            Key structure typically includes:

            -   `smiles` (str): The input target molecule SMILES.
            -   `steps` (list[dict]): A list of reaction steps. Each step includes:

                -   `step` (str): Step number.
                -   `reactants` (list[dict]): List of reactant molecules for the step.
                -   `products` (list[dict]): List of product molecules for the step.
                -   `reagents` (list[dict]): Predicted reagents for the step.
                -   `conditions` (dict): Predicted reaction conditions.
                -   `reactionmetrics` (list[dict]): Contains metrics like
                    `closestliterature` (predicted literature references).
                -   ... (other fields from `format_output` and `rec_run_prithvi`)

            -   `dependencies` (dict): Defines the relationship between steps.
            -   ... (potentially other top-level keys from `format_output`)

            Job-specific logs are created in `logs/<date>/job_<job_id>.log`.
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
            hallucination_check=hallucination_check)
        output_data = format_output(result_dict)
        output_data = add_metadata(output_data)
        return output_data
    finally:
        # Clean up handlers
        log._logger.removeHandler(handler)
        handler.close()
        context_logger.reset(token)


def add_metadata(output_data: dict) -> dict:
    """Adds predicted reagents, conditions, and literature to reaction steps.

    Iterates through each step in the `output_data` and calls the respective
    agents from `src.metadata` to predict and append:

    -   Reagents (to `step['reagents']`)
    -   Reaction conditions (to `step['conditions']`)
    -   Closest literature references (to `step['reactionmetrics'][0]['closestliterature']`)

    Args:
        output_data (dict):
            The formatted retrosynthesis output dictionary, typically from
            `format_output()`. It must contain a 'steps' key, which is a list
            of dictionaries. Each step dictionary in the list is expected to have:

            -   `reactants` (list[dict]): Input for `reagent_agent`, etc.
            -   `products` (list[dict]): Input for `reagent_agent`, etc.
            -   `reagents` (list): An existing list to which predicted reagents are appended.
            -   `conditions` (any): Will be overwritten by predicted conditions dict.
            -   `reactionmetrics` (list[dict]): A list where the first element is a
                dict that will have a `closestliterature` key added/updated.

    Returns:
        dict:
            The input `output_data` dictionary, modified in-place to include
            the predicted metadata within each step.
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
