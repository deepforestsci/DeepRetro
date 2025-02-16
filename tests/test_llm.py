import os
import ast

import rootutils
root_dir = rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

import pytest
from src.utils.llm import *
from src.variables import *

from tests.variables_test import *

def test_call_llm_SUCCESS():
    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING)
    assert status_code == 200

def test_call_llm_openai_SUCCESS():
    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")
    assert status_code == 200

def test_all_openai_SUCCESS():
    successful_tests = []
    failed_tests = []
    print("models: ", OPENAI_MODELS)
    for model in OPENAI_MODELS:
        try:
            status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING, LLM=model)
            assert status_code == 200
            successful_tests.append(model)
        except Exception as e:
            print(f"Error: {e}\n Model: {model}")
            failed_tests.append(model)
    if failed_tests:
        print(f"Failed models: {failed_tests}")

def test_call_llm_deepseek_SUCCESS():
    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")
    assert status_code == 200

def test_split_cot_json_SUCCESS():
    status_code, thinking_steps, json_content = split_cot_json(VALID_CLAUDE_RESPONSE)
    assert status_code == 200

def test_split_cot_json_FAIL_501():
    status_code, thinking_steps, json_content = split_cot_json(EMPTY_RESPONSE)
    assert status_code == 501

def test_split_json_openai_SUCCESS():
    status_code, res_OpenAI = call_LLM(molecule="CC(=O)CC", LLM="gpt-4o")
    assert status_code == 200

def test_split_json_openai_FAIL_502():
    status_code, res_OpenAI = split_json_openAI("")
    assert status_code == 502

def test_split_json_deepseek_SUCCESS():
    status_code, thinking_steps, json_content = split_json_deepseek(VALID_CLAUDE_RESPONSE)
    assert status_code == 200

def test_split_json_deepseek_FAIL_503():
    status_code, thinking_steps, json_content = split_json_deepseek(EMPTY_RESPONSE)
    assert status_code == 503


if __name__ == '__main__':
    pytest.main()
