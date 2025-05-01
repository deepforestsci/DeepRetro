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

# df = pd.read_csv(f"{root_dir}/data/challenge_mols.csv")
df = pd.read_csv(f"{root_dir}/data/USPTO_190_sample.csv")
column = "smiles"
if column not in df.columns:
    column = "input"
mols_dfs = df[column].to_list()
mapper = {
    'm0p0': 'claude-3-opus-20240229', 'm0p1': 'anthropic/claude-3-7-sonnet-20250219:adv',
    'm1p0': 'fireworks_ai/accounts/fireworks/models/deepseek-r1', 'm1p1': 'fireworks_ai/accounts/fireworks/models/deepseek-r1:adv',
}

folder_list = ["m0p1:Pistachio_100+", "m1p0:Pistachio_100+","m0p1:USPTO", "m1p0:USPTO", ]
MAX_CONCURRENT_TASKS = 4  # Process 3 molecules in parallel

async def process_molecule(mol, folder, run_no, mol_no):
    """Process a single molecule asynchronously."""
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
            
        output_path = f"{root_dir}/results/USPTO_190_sample/{folder}/run_{run_no}/{mol_no}.json"
        with open(output_path, "w") as f:
            json.dump(res_dict, f, indent=4)
            
    except Exception as e:
        print("Error in molecule:", mol)
        print("Error:", e)
        
    time2 = time.time()
    print(f"Time taken for {mol} on {folder}: {time2-time1} seconds")
    return mol

async def main_async():
    for run_no in range(0, 1):
        for folder in folder_list:
            # Create directories if they don't exist
            if not os.path.exists(f"{root_dir}/results/USPTO_190_sample/{folder}"):
                os.makedirs(f"{root_dir}/results/USPTO_190_sample/{folder}")
            if not os.path.exists(f"{root_dir}/results/USPTO_190_sample/{folder}/run_{run_no}"):
                os.makedirs(f"{root_dir}/results/USPTO_190_sample/{folder}/run_{run_no}")
            
            # Process molecules in batches to limit concurrency
            tasks = []
            for mol_no,mol in enumerate(mols_dfs):
                if len(tasks) >= MAX_CONCURRENT_TASKS:
                    # Wait for one task to complete when we reach the limit
                    done, pending = await asyncio.wait(tasks, return_when=asyncio.FIRST_COMPLETED)
                    tasks = list(pending)  # Keep only the pending tasks
                
                # Add a new task to the list
                task = asyncio.create_task(process_molecule(mol, folder, run_no, mol_no))
                tasks.append(task)
            
            # Wait for all remaining tasks to complete
            if tasks:
                await asyncio.gather(*tasks)

# Run the async main function
if __name__ == "__main__":
    asyncio.run(main_async())