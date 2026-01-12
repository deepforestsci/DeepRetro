# Single step retrosynthesis runner for USPTO-50k test dataset
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
from src.rec_prithvi import single_run_DeepRetro
from src.cache import clear_cache_for_molecule

# Load USPTO-50k test dataset (250 molecules subset)
df = pd.read_csv(f"{root_dir}/data/uspto_50k_test_250.csv")
mols_dfs = df['input'].to_list()  # Using 'input' column for target molecules

# Model mapper
mapper_new = {
    'claude3opus': 'claude-3-opus-20240229:adv',
    'deepseekr1': 'fireworks_ai/accounts/fireworks/models/deepseek-r1',
    'claude35sonnet': 'anthropic/claude-3-5-sonnet-20241022:adv',
    'claude37sonnet': 'anthropic/claude-3-7-sonnet-20250219:adv',
    'claude4sonnet': 'anthropic/claude-sonnet-4-20250514:adv',
    'claude4opus': 'anthropic/claude-opus-4-20250514:adv',
    'claude45sonnet': 'anthropic/claude-sonnet-4-5-20250929:adv',
}

folder_list = [
    "claude3opus:Pistachio_100+",
    # "deepseekr1:Pistachio_100+",
    "claude4sonnet:Pistachio_100+",
    "claude4opus:Pistachio_100+",
    # "claude35sonnet:Pistachio_100+",
    "claude37sonnet:Pistachio_100+",
    "claude45sonnet:Pistachio_100+",
]
MAX_CONCURRENT_TASKS = 5  # Process 5 molecules in parallel

results_dir = f"{root_dir}/results/single_step"


async def process_molecule(mol, folder, run_no):
    """Process a single molecule asynchronously using single_run_DeepRetro."""
    molecule = mol
    llm = mapper_new[folder.split(":")[0]]
    az_model = folder.split(":")[1]
    time1 = time.time()

    try:
        clear_cache_for_molecule(molecule)
        print(f"Running {molecule} with {llm} and {az_model}")

        # Run the potentially blocking single_run_DeepRetro function in a thread pool
        loop = asyncio.get_event_loop()
        with ThreadPoolExecutor() as pool:
            res_dict, solved = await loop.run_in_executor(
                pool,
                functools.partial(single_run_DeepRetro,
                                  molecule,
                                  llm=llm,
                                  az_model=az_model,
                                  hallucination_check="True",
                                  stability_flag="True"))

        # Add solved status to the result dictionary
        result_output = {
            "molecule": molecule,
            "solved": solved,
            "result": res_dict
        }

        output_path = f"{results_dir}/{folder}/run_{run_no}/{mol}.json"
        with open(output_path, "w") as f:
            json.dump(result_output, f, indent=4)

    except Exception as e:
        print("Error in molecule:", mol)
        print("Error:", e)

    time2 = time.time()
    print(f"Time taken for {mol} on {folder}: {time2-time1} seconds")
    return mol


async def main_async():
    for run_no in range(1, 2):
        for folder in folder_list:
            # Create directories if they don't exist
            if not os.path.exists(f"{results_dir}/{folder}"):
                os.makedirs(f"{results_dir}/{folder}")
            if not os.path.exists(f"{results_dir}/{folder}/run_{run_no}"):
                os.makedirs(f"{results_dir}/{folder}/run_{run_no}")

            # Process molecules in batches to limit concurrency
            tasks = []
            for mol in mols_dfs:
                if len(tasks) >= MAX_CONCURRENT_TASKS:
                    # Wait for one task to complete when we reach the limit
                    done, pending = await asyncio.wait(
                        tasks, return_when=asyncio.FIRST_COMPLETED)
                    tasks = list(pending)  # Keep only the pending tasks

                # Add a new task to the list
                task = asyncio.create_task(
                    process_molecule(mol, folder, run_no))
                tasks.append(task)

            # Wait for all remaining tasks to complete
            if tasks:
                await asyncio.gather(*tasks)


# Run the async main function
if __name__ == "__main__":
    asyncio.run(main_async())

