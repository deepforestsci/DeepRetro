import os
import ast
import pytest
from unittest.mock import patch

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)


from src.utils.llm import call_LLM
from tests.variables_test import VALID_SMILE_STRING, DEEPSEEK_FIREWORKS_MODEL, CLAUDE_ADV_MODEL


@patch("src.utils.llm.call_LLM")
def test_claude_adv_success(mock_call_llm):
    mock_call_llm.return_value = (200, "mocked response")
    status, response = call_LLM(molecule=VALID_SMILE_STRING, LLM=CLAUDE_ADV_MODEL)
    assert status == 200

# OpenAI tests are commented, because OpenAI models are not being used.
# def test_openai_adv_success():

#     status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING,
#                                      LLM=OPENAI_MODEL_ADV)
#     if not res_text:
#         assert status_code == 404
#     assert status_code == 200


@patch("src.utils.llm.call_LLM")
def test_deepseek_adv_success(mock_call_llm):
    mock_call_llm.return_value = (200, "mocked response")
    status, response = call_LLM(molecule="CC1=NC=C(N1CCO)[N+]([O-])=O", LLM=DEEPSEEK_FIREWORKS_MODEL)
    assert status == 200


if __name__ == '__main__':
    pytest.main()
