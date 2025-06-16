"""Processes a sample of 500 entries from the USPTO-50k test dataset.

This script first loads the `test.csv` file from the USPTO-50k dataset,
(Note: The path to this CSV is currently hardcoded to
`/Users/shreyasv/Desktop/Deepchem/data/USPTO-50k/test.csv` and may need adjustment
for different environments).

It then takes a random sample of 500 entries from this dataset, sorts them by
their original index, and processes each entry.
For each sampled entry, the `run_prithvi` function (imported from a local `utils`
module) is called with the input from the row. The resulting dictionary is saved
as a JSON file in the `results/USPTO-50k_500/` directory, named by the
original index of the row.

Errors during the processing of any entry are caught and printed to standard output.
The script also imports `clear_cache` from `utils` but does not appear to use it.

Dependencies:
    - pandas: For reading and processing the CSV file.
    - A local `utils.py` module containing `run_prithvi` and `clear_cache`.

Assumed File Structure:
    - The USPTO-50k `test.csv` at the specified hardcoded path.
    - A directory `results/USPTO-50k_500/` for saving output JSON files.
"""
# load up USPTO-50k test dataset
import pandas as pd
from utils import clear_cache, run_prithvi
import json

df = pd.read_csv("/Users/shreyasv/Desktop/Deepchem/data/USPTO-50k/test.csv")

# now extract 500 entries from the dataset
df = df.sample(500, random_state=42)
df.sort_index(inplace=True)
# df.to_csv("USPTO-50k_500.csv")

# now run prithvi on the dataset

for i, row in df.iterrows():
    print(f"Running {i}")
    try:
        res_dict = run_prithvi(row['input'])
        with open(f"results/USPTO-50k_500/{i}.json", "w") as f:
            json.dump(res_dict, f, indent=4)
    except Exception as e:
        print("Error in molecule:", i)
        print("Error:", e)
