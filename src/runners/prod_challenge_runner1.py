"""Runs retrosynthesis predictions for a subset of failed cases from USPTO-50k.

This script reads a predefined CSV file (`USPTO-50k_500_failed.csv`), specifically
the first 15 rows after reading. For each row, it extracts an input (presumably
a SMILES string or similar molecular representation) and processes it using the
`run_prithvi` function (imported from a `utils` module, expected to be in the
same directory or Python path).

The results from `run_prithvi` for each input are saved as individual JSON files
in the `results/USPTO-50k_500/` directory, named by the index of the row from
the input DataFrame. Errors during processing are caught and printed to standard
output.

Dependencies:
    - pandas: For reading and processing the CSV file.
    - tqdm: For displaying a progress bar during processing.
    - A local `utils.py` module containing the `run_prithvi` function.

Assumed File Structure:
    - `USPTO-50k_500_failed.csv` in the same directory as the script.
    - A directory `results/USPTO-50k_500/` for saving output JSON files.
"""
import pandas as pd
from utils import run_prithvi
import json
from tqdm import tqdm

df = pd.read_csv("USPTO-50k_500_failed.csv", index_col=0)
df = df[:15]


for i, row in tqdm(df.iterrows(), total=len(df)):
    # print(f"Running {i}")
    try:
        res_dict = run_prithvi(row['input'])
        with open(f"results/USPTO-50k_500/{i}.json", "w") as f:
            json.dump(res_dict, f, indent=4)
    except Exception as e:
        print("Error in molecule:", i)
        print("Error:", e)