import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

from dotenv import load_dotenv
load_dotenv()

from src.utils.llm import call_LLM, split_cot_json, split_json_deepseek

SMILE_STRING = "CC(=O)O"


def test_call_llm_success():

    from tests.variables_test import VALID_SMILE_STRING

    status_code, _ = call_LLM(molecule=VALID_SMILE_STRING)
    assert status_code == 200


def test_split_cot_json_success():
    
    from tests.variables_test import VALID_CLAUDE_RESPONSE

    status_code, _, _ = split_cot_json(VALID_CLAUDE_RESPONSE)
    assert status_code == 200


def test_split_cot_json_fail_501():
    
    from tests.variables_test import EMPTY_RESPONSE
    
    status_code, _, _ = split_cot_json(EMPTY_RESPONSE)
    assert status_code == 501


def test_call_llm_deepseek_success():
    '''
    Testing deepseek model, takes too long to run
    '''
    from tests.variables_test import VALID_SMILE_STRING, DEEPSEEK_FIREWORKS_MODEL

    status_code, res_text = call_LLM(molecule="CC(=O)O", LLM=DEEPSEEK_FIREWORKS_MODEL)
    print(res_text)
    assert status_code == 200


def test_split_json_deepseek_success():

    from tests.variables_test import DEEPSEEK_ADV_VALID_RESPONSE
    status_code, _, _ = split_json_deepseek(DEEPSEEK_ADV_VALID_RESPONSE)
    assert status_code == 200


def test_split_json_deepseek_fail_503():
    from tests.variables_test import EMPTY_RESPONSE
    
    status_code, _, _ = split_json_deepseek(EMPTY_RESPONSE)
    assert status_code == 503


def test_split_json_master_success():
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL
    from src.utils.llm import split_json_master

    status_code, _, json_content = split_json_master(VALID_CLAUDE_RESPONSE, model=CLAUDE_MODEL)
    assert status_code == 200


def test_split_json_master_fail():
    from tests.variables_test import EMPTY_RESPONSE, CLAUDE_MODEL
    from src.utils.llm import split_json_master

    status_code, _, json_content = split_json_master(EMPTY_RESPONSE, model=CLAUDE_MODEL)
    assert status_code != 200


def test_validate_json_success_200():
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL
    from src.utils.llm import validate_split_json, split_json_master

    status_code, _, json_content = split_json_master(
        VALID_CLAUDE_RESPONSE, CLAUDE_MODEL)

    status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    assert status_code == 200


def test_validate_json_fail():
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL, EMPTY_RESPONSE
    from src.utils.llm import validate_split_json, split_json_master

    status_code, _, json_content = split_json_master(
        EMPTY_RESPONSE, CLAUDE_MODEL)

    status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    assert status_code != 200

def test_validity_check_success():
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL, VALID_SMILE_STRING
    from src.utils.llm import validate_split_json, split_json_master, validity_check

    status_code, _, json_content = split_json_master(
        VALID_CLAUDE_RESPONSE, CLAUDE_MODEL)

    status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    output_pathways, output_explanations, output_confidence = validity_check(
        VALID_SMILE_STRING, res_molecules, res_explanations, res_confidence)

    assert output_pathways

def test_validity_check_fail():
    from tests.variables_test import CLAUDE_MODEL, EMPTY_RESPONSE
    from src.utils.llm import validate_split_json, split_json_master, validity_check

    status_code, _, json_content = split_json_master(
        EMPTY_RESPONSE, CLAUDE_MODEL)

    status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    output_pathways, output_explanations, output_confidence = validity_check(
        "", res_molecules, res_explanations, res_confidence)

    assert not output_pathways


# OpenAI tests

# def test_call_llm_openai_success():

#     from tests.variables_test import VALID_SMILE_STRING

#     status_code, _ = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")
#     assert status_code == 200

# def test_all_openai_success():

#     from tests.variables_test import VALID_SMILE_STRING

#     successful_tests = []
#     failed_tests = []
#     print("models: ", OPENAI_MODELS)
#     for model in OPENAI_MODELS:
#         try:
#             status_code, res_text = call_LLM(molecule=VALID_SMILE_STRING, LLM=model)
#             assert status_code == 200
#             successful_tests.append(model)
#         except Exception as e:
#             print(f"Error: {e}\n Model: {model}")
#             failed_tests.append(model)
#     if failed_tests:
#         print(f"Failed models: {failed_tests}")

# def test_split_json_openai_success():
    
#     from tests.variables_test import VALID_SMILE_STRING

#     status_code, _ = call_LLM(molecule=VALID_SMILE_STRING, LLM="gpt-4o")
#     assert status_code == 200

# def test_split_json_openai_fail_502():

#     status_code, _ = split_json_openAI("")
#     assert status_code == 502

if __name__ == '__main__':
    pytest.main()
