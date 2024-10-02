import pandas as pd
from utils import run_prithvi
import json
from tqdm import tqdm

df = pd.read_csv("USPTO-50k_500_failed.csv", index_col=0)
df = df[60:]


for i, row in tqdm(df.iterrows(), total=len(df)):
    # print(f"Running {i}")
    try:
        res_dict = run_prithvi(row['input'])
        with open(f"results/USPTO-50k_500/{i}.json", "w") as f:
            json.dump(res_dict, f, indent=4)
    except Exception as e:
        print("Error in molecule:", i)
        print("Error:", e)