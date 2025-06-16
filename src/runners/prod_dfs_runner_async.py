"""Asynchronously runs a series of Depth-First Search (DFS) experiments.

This script is an asynchronous version of `prod_dfs_runner.py`. It uses `asyncio`
and `concurrent.futures.ThreadPoolExecutor` to process multiple molecules in
parallel, aiming to speed up the overall experimental runs.

The core workflow, including loading molecules, mapping configurations, defining
experimental folders, and iterating through run numbers, is similar to the
synchronous version. The main difference lies in how individual molecule processing
is handled.

Key Features:
- Asynchronous processing of molecules using `asyncio`.
- Limits the number of concurrently processed molecules via `MAX_CONCURRENT_TASKS`.
- The synchronous `src.main.main()` function is run in a `ThreadPoolExecutor` to
  avoid blocking the asyncio event loop.
- Saves results to JSON files with a `_hallucination.json` suffix in a structured
  directory: `results/dfs/<condition_folder>/run_<run_no>/`.

Dependencies:
    - (Same as `prod_dfs_runner.py` plus `asyncio`, `concurrent.futures`).
"""
# load up USPTO-50k test dataset
import rootutils
import os
import pandas as pd
import json
import time
import asyncio
import functools
from concurrent.futures import ThreadPoolExecutor

root_dir = rootutils.setup_root(".",
                               indicator=".project-root",
                               pythonpath=True)
from src.main import main
from src.cache import clear_cache_for_molecule

df = pd.read_csv(f"{root_dir}/results/dfs/dataset.csv")
mols_dfs = df['smiles'].to_list()
mapper = {
    'm0p0': 'claude-3-opus-20240229', 'm0p1': 'anthropic/claude-3-7-sonnet-20250219:adv',
    'm1p0': 'fireworks_ai/accounts/fireworks/models/deepseek-r1', 'm1p1': 'fireworks_ai/accounts/fireworks/models/deepseek-r1:adv',
}

folder_list = ["m0p1:Pistachio_100+", "m1p0:Pistachio_100+"]
MAX_CONCURRENT_TASKS = 6  # Process 5 molecules in parallel

async def process_molecule(mol, folder, run_no):
    """Asynchronously processes a single molecule for a given experimental setup.

    This function handles the logic for one molecule: clearing its cache, determining
    LLM and AZ model settings, running the main processing logic via `src.main.main`
    (offloaded to a thread pool), and saving the results to a JSON file.

    Args:
        mol (str): The SMILES string of the molecule to process.
        folder (str): The string representing the current experimental condition
            (e.g., "m0p1:Pistachio_100+").
        run_no (int): The current run number (for multiple trials).

    Returns:
        str: The SMILES string of the processed molecule (primarily for task tracking).
    """
    molecule = mol
    llm = mapper[folder.split(":")[0]]
    az_model = folder.split(":")[1]
    time1 = time.time()
    
    try:
        clear_cache_for_molecule(molecule)
        print(f"Running {molecule} with {llm} and {az_model}")
        
        # Run the potentially blocking main function in a thread pool
        loop = asyncio.get_event_loop()
        with ThreadPoolExecutor() as pool:
            res_dict = await loop.run_in_executor(
                pool, 
                functools.partial(main, molecule, llm=llm, az_model=az_model, hallucination_check="True", stability_flag="True")
            )
            
        output_path = f"{root_dir}/results/dfs/{folder}/run_{run_no}/{mol}_hallucination.json"
        with open(output_path, "w") as f:
            json.dump(res_dict, f, indent=4)
            
    except Exception as e:
        print("Error in molecule:", mol)
        print("Error:", e)
        
    time2 = time.time()
    print(f"Time taken for {mol} on {folder}: {time2-time1} seconds")
    return mol

async def main_async():
    """Main asynchronous function to orchestrate experimental runs.

    Iterates through run numbers and experimental folder configurations. For each
    combination, it creates and manages asynchronous tasks for processing all
    molecules in the dataset, respecting the `MAX_CONCURRENT_TASKS` limit.
    Ensures that necessary output directories are created.
    """
    for run_no in range(9, 12):
        for folder in folder_list:
            # Create directories if they don't exist
            if not os.path.exists(f"{root_dir}/results/dfs/{folder}"):
                os.makedirs(f"{root_dir}/results/dfs/{folder}")
            if not os.path.exists(f"{root_dir}/results/dfs/{folder}/run_{run_no}"):
                os.makedirs(f"{root_dir}/results/dfs/{folder}/run_{run_no}")
            
            # Process molecules in batches to limit concurrency
            tasks = []
            for mol in mols_dfs:
                if len(tasks) >= MAX_CONCURRENT_TASKS:
                    # Wait for one task to complete when we reach the limit
                    done, pending = await asyncio.wait(tasks, return_when=asyncio.FIRST_COMPLETED)
                    tasks = list(pending)  # Keep only the pending tasks
                
                # Add a new task to the list
                task = asyncio.create_task(process_molecule(mol, folder, run_no))
                tasks.append(task)
            
            # Wait for all remaining tasks to complete
            if tasks:
                await asyncio.gather(*tasks)

# Run the async main function
if __name__ == "__main__":
    asyncio.run(main_async())