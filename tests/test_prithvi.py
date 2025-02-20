import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from src.utils.llm import call_LLM, split_cot_json, split_json_openAI, split_json_deepseek
from src.variables import OPENAI_MODELS
from tests.variables_test import VALID_SMILE_STRING

DEEPSEEK_MODEL_ADV = "deepinfra/deepseek-ai/DeepSeek-R1:adv"
OPENAI_MODEL_ADV = "claude-3-opus-20240229:adv"
CLAUDE_MODEL_ADV = "claude-3-opus-20240229:adv"

print("Used smile string: ", VALID_SMILE_STRING)






if __name__ == '__main__':
    pytest.main()
