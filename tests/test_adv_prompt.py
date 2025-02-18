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


def test_claude_adv_success():

    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING,
                                     LLM="claude-3-opus-20240229:adv")
    if not res_text:
        assert status_code == 404
    assert status_code == 200



if __name__ == '__main__':
    pytest.main()