"""Runs a series of Depth-First Search (DFS) experiments with varied configurations.

This script is designed to automate running the main processing pipeline
(via `src.main.main`) for a list of molecules under different experimental
conditions. The conditions are defined by combinations of Large Language Models (LLMs)
and AiZynthFinder (AZ) model settings.

Workflow:
1.  Loads a list of molecules (SMILES strings) from `results/dfs/dataset.csv`.
2.  Defines a `mapper` to translate short configuration codes (e.g., "m0p0") to
    specific LLM model identifiers.
3.  Defines a `folder_list` where each item represents an experimental condition,
    combining an LLM configuration code and an AZ model/dataset identifier
    (e.g., "m0p1:Pistachio_100+").
4.  Iterates through a specified range of run numbers (e.g., for multiple trials).
5.  For each run number and each experimental condition in `folder_list`:
    a.  Ensures necessary output directories are created under `results/dfs/`.
    b.  For every molecule in the loaded dataset:
        i.  Clears any existing cache entries for that specific molecule using
            `clear_cache_for_molecule`.
        ii. Extracts the LLM and AZ model settings from the current folder/condition.
        iii. Calls `src.main.main()` with the molecule, LLM, AZ model, and with
             `stability_flag` and `hallucination_check` set to "True".
        iv. Saves the output dictionary as a JSON file in a structured directory:
            `results/dfs/<condition_folder>/run_<run_no>/<molecule_smiles>_stability_hallucination.json`.
        v.  Prints the time taken for processing each molecule under the current condition.

This script facilitates batch processing and systematic evaluation of the main
retrosynthesis pipeline across different models and settings.

Dependencies:
    - rootutils: For project root setup.
    - pandas: For reading the input dataset CSV.
    - src.main: Provides the `main()` processing function.
    - src.cache: Provides `clear_cache_for_molecule()`.
    - os, json, time: Standard Python libraries.
"""
# load up USPTO-50k test dataset
import rootutils
import os
import pandas as pd
import json
import time
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

# # now run prithvi on the dataset
# folder_list = [ "m1p0:Pistachio_100+", "m1p0:Pistachio_50", "m0p1:Pistachio_100+", "m0p1:Pistachio_50"]
folder_list = [ "m0p1:Pistachio_100+", "m1p0:Pistachio_100+"]

for run_no in range(9,12):
    for folder in folder_list:
        if not os.path.exists(f"{root_dir}/results/dfs/{folder}"):
            os.makedirs(f"{root_dir}/results/dfs/{folder}")
        if not os.path.exists(f"{root_dir}/results/dfs/{folder}/run_{run_no}"):
            os.makedirs(f"{root_dir}/results/dfs/{folder}/run_{run_no}")
        for mol in mols_dfs:
            molecule = mol
            print(f"Running {mol}")
            llm = mapper[folder.split(":")[0]]
            az_model = folder.split(":")[1]
            time1 = time.time()
            try:
                clear_cache_for_molecule(molecule)
                print(f"Running {molecule} with {llm} and {az_model}")
                res_dict = main(molecule ,llm=llm, az_model=az_model, stability_flag="True", hallucination_check="True")
                with open(f"{root_dir}/results/dfs/{folder}/run_{run_no}/{mol}_stability_hallucination.json", "w") as f:
                    json.dump(res_dict, f, indent=4)
            except Exception as e:
                print("Error in molecule:", mol)
                print("Error:", e)
            time2 = time.time()
            print(f"Time taken for {mol} on {folder}: {time2-time1} seconds")
