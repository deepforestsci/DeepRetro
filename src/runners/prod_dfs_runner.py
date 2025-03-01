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
