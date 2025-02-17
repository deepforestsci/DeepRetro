# load up USPTO-50k test dataset
import rootutils
import os
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)
from src.main import main
from src.cache import clear_cache_for_molecule
import pandas as pd
import json

df = pd.read_csv(f"{root_dir}/results/dfs/dataset.csv")
mols_dfs = df['smiles'].to_list()

# # now run prithvi on the dataset
folder = "m1p1"
if not os.path.exists(f"{root_dir}/results/dfs/{folder}"):
    os.makedirs(f"{root_dir}/results/dfs/{folder}")
for mol in mols_dfs:
    molecule = mol
    print(f"Running {mol}")
    # try:
    #     res_dict = run_prithvi(molecule)
    #     with open(f"{root_dir}/results/mols_hard/{mol}.json", "w") as f:
    #         json.dump(res_dict, f, indent=4)
    # except Exception as e:
    #     print("Error in molecule:", mol)
    #     print("Error:", e)
    try:
        clear_cache_for_molecule(molecule)
        # "azure_ai/DeepSeek-R1", "deepinfra/deepseek-ai/DeepSeek-R1"
        res_dict = main(molecule ,llm="deepinfra/deepseek-ai/DeepSeek-R1:adv")
        with open(f"{root_dir}/results/dfs/{folder}/{mol}.json", "w") as f:
            json.dump(res_dict, f, indent=4)
    except Exception as e:
        print("Error in molecule:", mol)
        print("Error:", e)
