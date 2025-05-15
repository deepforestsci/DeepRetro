import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from dotenv import load_dotenv
load_dotenv()

from src.utils.llm import call_LLM, split_cot_json, split_json_deepseek

SMALL_SMILE_STRING = "CC(=O)O"
LARGE_SMILE_STRING = "CC(N)C(=O)NC1=C(C)C=CC=C1C"

def test_call_llm_success():
    """Tests call_LLM function with valid smile string.
    
    Expected output:
        status_code: 200.
        res_text: str.
    """
    status_code, res_text = call_LLM(molecule=LARGE_SMILE_STRING)
    
    assert status_code == 200
    assert isinstance(res_text, str)


def test_split_cot_json_success():
    """Tests split_cot_json function with valid response. split_cot_json() splits the response based on <cot>, <thinking> and <json> tag.
    
    Expected output:
        status_code: 200.
        thinking_steps: list containing items
        json_content: str.
    """
    from tests.variables_test import VALID_CLAUDE_RESPONSE

    status_code, thinking_steps, json_content = split_cot_json(VALID_CLAUDE_RESPONSE)
    
    assert status_code == 200
    assert isinstance(thinking_steps, list)
    assert thinking_steps
    assert isinstance(json_content, str)
    assert json_content


def test_split_cot_json_fail_501():
    """Tests split_cot_json function with empty response.
    
    Expected output:
        status_code: 501
        thinking_steps: []
        json_content: ""
    """
    from tests.variables_test import EMPTY_RESPONSE
    
    status_code, thinking_steps, json_content = split_cot_json(EMPTY_RESPONSE)
    
    assert status_code == 501
    assert thinking_steps == []
    assert json_content == ""


def test_call_llm_deepseek_success():
    '''Testing deepseek model, hosted in fireworks.
    
    Expected output:
        status_code: 200.
        res_text: str.
    '''
    from tests.variables_test import DEEPSEEK_FIREWORKS_MODEL

    status_code, res_text = call_LLM(molecule=LARGE_SMILE_STRING, LLM=DEEPSEEK_FIREWORKS_MODEL)
    
    assert status_code == 200
    assert isinstance(res_text, str)
    assert res_text


def test_split_json_deepseek_success():
    """Tests split_json_deepseek function with valid response. split_json_deepseek splits the response based on <thinking> and <jon> tag.
    
    Expected output:
        status_code: 200.
        thinking_steps: list containing items
        json_content: str
    """
    from tests.variables_test import DEEPSEEK_ADV_VALID_RESPONSE
    
    status_code, thinking_steps, json_content = split_json_deepseek(DEEPSEEK_ADV_VALID_RESPONSE)
    
    assert status_code == 200
    assert isinstance(thinking_steps, list)
    assert isinstance(json_content, str)


def test_split_json_deepseek_fail_503():
    """Tests split_json_deepseek function with empty response.
    
    Expected output:
        status_code: 503
        thinking_step: []
        json_content: ""
    """
    from tests.variables_test import EMPTY_RESPONSE
    
    status_code, thinking_steps, json_content = split_json_deepseek(EMPTY_RESPONSE)
    
    assert status_code == 503
    assert thinking_steps == []
    assert json_content == ""


if __name__ == '__main__':
    pytest.main()
