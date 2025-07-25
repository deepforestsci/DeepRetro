import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)


from src.utils.llm import call_LLM
from tests.variables_test import VALID_SMILE_STRING, DEEPSEEK_FIREWORKS_MODEL, CLAUDE_ADV_MODEL


def test_claude_adv_success():
    """Tests call_LLM function with advance claude model.
    """
    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING,
                                     LLM=CLAUDE_ADV_MODEL)
    if not res_text:
        print("res_text is empty")
        status_code = 400
    
    assert status_code == 200

# OpenAI tests are commented, because OpenAI models are not being used.
# def test_openai_adv_success():

#     status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING,
#                                      LLM=OPENAI_MODEL_ADV)
#     if not res_text:
#         assert status_code == 404
#     assert status_code == 200


def test_deepseek_adv_success():
    """Tests call_LLM function with advance deepseek model(hosted in fireworks).
    """

    status_code, res_text = call_LLM(molecule="CC1=NC=C(N1CCO)[N+]([O-])=O",
                                     LLM=DEEPSEEK_FIREWORKS_MODEL)
    if not res_text:
        print("res_text is empty")
        status_code = 400

    assert status_code == 200



if __name__ == '__main__':
    pytest.main()
