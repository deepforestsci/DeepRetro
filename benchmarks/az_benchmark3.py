import rootutils

root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tqdm import tqdm
import os
import json
from src.utils.parse import format_output
from src.utils.az import run_az

import pandas as pd

df = pd.read_csv("/home/ubuntu/recursiveLLM/benchmarks/sampled_molecules.csv")

az_model_list = ["Pistachio_100+", "USPTO"]
for az_model in az_model_list:
    if not os.path.exists(
            f"/home/ubuntu/recursiveLLM/results/benchmarking/AZ_paper/{az_model}"
    ):
        os.makedirs(
            f"/home/ubuntu/recursiveLLM/results/benchmarking/AZ_paper/{az_model}"
        )
    print(f"Running {az_model}")
    for i, row in tqdm(df.iterrows(), total=len(df)):
        if i >= 126 and i < 189:
            print(f"Running {i}")
            solved, result_dict = run_az(smiles=row['input'], az_model=az_model)
            output_data = format_output(result_dict[0])
            try:
                with open(
                        f"/home/ubuntu/recursiveLLM/results/benchmarking/AZ_paper/{az_model}/{i}.json",
                        "w") as f:
                    json.dump(output_data, f, indent=4)
            except Exception as e:
                print("Error in molecule:", i)
                print("Error:", e)
