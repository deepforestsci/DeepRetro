import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

from src.utils.llm import call_LLM, split_cot_json, split_json_openAI, split_json_deepseek
from src.variables import OPENAI_MODELS

    

def test_call_llm_success():

    from tests.variables_test import VALID_SMILE_STRING

    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING)
    assert status_code == 200

def test_call_llm_openai_success():

    from tests.variables_test import VALID_SMILE_STRING

    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")
    assert status_code == 200

def test_all_openai_success():

    from tests.variables_test import VALID_SMILE_STRING

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

def test_call_llm_deepseek_success():

    from tests.variables_test import VALID_SMILE_STRING

    status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")
    assert status_code == 200

def test_split_cot_json_success():
    
    from tests.variables_test import VALID_CLAUDE_RESPONSE

    status_code, thinking_steps, json_content = split_cot_json(VALID_CLAUDE_RESPONSE)
    assert status_code == 200

def test_split_cot_json_fail_501():
    
    from tests.variables_test import EMPTY_RESPONSE
    
    status_code, thinking_steps, json_content = split_cot_json(EMPTY_RESPONSE)
    assert status_code == 501

def test_split_json_openai_success():
    
    from tests.variables_test import VALID_SMILE_STRING

    status_code, _ = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")

    assert status_code == 200

def test_split_json_openai_fail_502():
    status_code, _ = split_json_openAI("")
    assert status_code == 502

def test_split_json_deepseek_success():

    from tests.variables_test import VALID_CLAUDE_RESPONSE
    status_code, thinking_steps, json_content = split_json_deepseek(VALID_CLAUDE_RESPONSE)
    assert status_code == 200

def test_split_json_deepseek_fail_503():
    from tests.variables_test import EMPTY_RESPONSE
    
    status_code, thinking_steps, json_content = split_json_deepseek(EMPTY_RESPONSE)
    assert status_code == 503


if __name__ == '__main__':
    pytest.main()
