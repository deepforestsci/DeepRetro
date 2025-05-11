import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root("..", indicator=".project-root", pythonpath=True)

from dotenv import load_dotenv
load_dotenv()


SMALL_SMILE_STRING = "CC(=O)O"
LARGE_SMILE_STRING = "CC(N)C(=O)NC1=C(C)C=CC=C1C"


def test_split_json_master_success():
    """Tests split_json_master function with valid response. split_json_master call the split_cot_json, split_json_deepseek and split_json_openAI functions, based on the model passed as argument.
    
    Expected output:
        status_code: 200.
        thinking_steps: list, containing items.
        json_content: str
    """
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL
    from src.utils.llm import split_json_master

    status_code, thinking_steps, json_content = split_json_master(VALID_CLAUDE_RESPONSE, model=CLAUDE_MODEL)
    
    assert status_code == 200
    assert isinstance(thinking_steps, list)
    assert thinking_steps
    assert isinstance(json_content, str)
    assert json_content


def test_split_json_master_fail():
    """Tests split_json_master function with empty response.

    Expected output:
        status_code: 501
        thinking_step: []
        json_content: ""
    """
    from tests.variables_test import EMPTY_RESPONSE, CLAUDE_ADV_MODEL
    from src.utils.llm import split_json_master

    status_code, thinking_steps, json_content = split_json_master(EMPTY_RESPONSE, model=CLAUDE_ADV_MODEL)
    assert status_code == 501
    assert thinking_steps == []
    assert json_content == ""


def test_validate_json_success_200():
    """Tests validate_split_json function with valid response from call_LLM. validate_split_json() extracts res_molecules, res_explanations, res_confidence from the json_content of the response.

    Expected output:
        status_code: 200
        res_molecules: list containing items.
        res_explanations: list containing items.
        res_confidence: list containing items.
    """
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_ADV_MODEL
    from src.utils.llm import validate_split_json, split_json_master

    status_code, _, json_content = split_json_master(
        VALID_CLAUDE_RESPONSE, CLAUDE_ADV_MODEL)

    status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    assert status_code == 200
    assert isinstance(res_molecules, list)
    assert res_molecules
    assert isinstance(res_explanations, list)
    assert res_explanations
    assert isinstance(res_confidence, list)
    assert res_confidence


def test_validate_json_fail():
    """Tests validate_split_json with empty response.

    Expected Output:
        status_code: 504
        res_molecules: []
        res_explanations: []
        res_confidence: []
    """
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL, EMPTY_RESPONSE
    from src.utils.llm import validate_split_json, split_json_master

    status_code, _, json_content = split_json_master(
        EMPTY_RESPONSE, CLAUDE_MODEL)

    status_code, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    assert status_code == 504
    assert res_molecules == []
    assert res_explanations == []
    assert res_confidence == []


def test_validity_check_success():
    """Tests validity_check with valid smile string. validity_check checks the validity of the molecules obtained from LLM.

    Expected Output:
        output_pathways: list containing output pathways.
        output_explanations: list containing explanations to output pathways.
        output_confidence: list containing confidence score.
    """
    from tests.variables_test import VALID_CLAUDE_RESPONSE, CLAUDE_MODEL, VALID_SMILE_STRING
    from src.utils.llm import validate_split_json, split_json_master, validity_check

    _, _, json_content = split_json_master(
        VALID_CLAUDE_RESPONSE, CLAUDE_MODEL)

    _, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    output_pathways, output_explanations, output_confidence = validity_check(
        VALID_SMILE_STRING, res_molecules, res_explanations, res_confidence)

    assert isinstance(output_pathways, list)
    assert output_pathways
    assert isinstance(output_explanations, list)
    assert output_explanations
    assert isinstance(output_confidence, list)
    assert output_confidence


def test_validity_check_fail():
    """Tests validity_check with empty response.

    Expected Outout:
        output_pathways: []
        output_explanations: []
        output_confidence: []
    """
    from tests.variables_test import CLAUDE_MODEL, EMPTY_RESPONSE
    from src.utils.llm import validate_split_json, split_json_master, validity_check

    _, _, json_content = split_json_master(
        EMPTY_RESPONSE, CLAUDE_MODEL)

    _, res_molecules, res_explanations, res_confidence = validate_split_json(
        json_content)

    output_pathways, output_explanations, output_confidence = validity_check(
        "", res_molecules, res_explanations, res_confidence)

    assert output_pathways == []
    assert output_explanations == []
    assert output_confidence == []


if __name__ == '__main__':
    pytest.main()
