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
