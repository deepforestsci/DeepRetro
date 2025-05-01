from litellm import completion
from typing import Optional
import litellm
from dotenv import load_dotenv
import rootutils

rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

from src.variables import SYS_PROMPT, USER_PROMPT
from src.prithvi import add_metadata

file_path = "/home/ubuntu/recursiveLLM/results/SICIC/mol3/mol3_path_6.json"

# read the json content
import json

with open(file_path) as f:
    output_data = json.load(f)

new_output_data = add_metadata(output_data)

# save the new output data
file_path_new = "/home/ubuntu/recursiveLLM/results/SICIC/mol3/mol3_path_6_new.json"
with open(file_path_new, 'w') as f:
    json.dump(new_output_data, f, indent=4)