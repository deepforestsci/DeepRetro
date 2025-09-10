import os
import ast
import pytest

import rootutils

root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from dotenv import load_dotenv

load_dotenv()

from src.utils.llm import call_LLM, split_cot_json, split_json_deepseek
from unittest.mock import patch, MagicMock

SMALL_SMILE_STRING = "CC(=O)O"
LARGE_SMILE_STRING = "CC(N)C(=O)NC1=C(C)C=CC=C1C"


def test_call_llm_success():
    """Tests call_LLM function with valid smile string.
    
    Expected output:
        status_code: 200.
        res_text: str.
    """
    from tests.variables_test import VALID_SMILE_STRING

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

    status_code, thinking_steps, json_content = split_cot_json(
        VALID_CLAUDE_RESPONSE)

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

    status_code, res_text = call_LLM(molecule=LARGE_SMILE_STRING,
                                     LLM=DEEPSEEK_FIREWORKS_MODEL)

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

    status_code, thinking_steps, json_content = split_json_deepseek(
        DEEPSEEK_ADV_VALID_RESPONSE)

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

    status_code, thinking_steps, json_content = split_json_deepseek(
        EMPTY_RESPONSE)

    assert status_code == 503
    assert thinking_steps == []
    assert json_content == ""


# OpenAI tests
# Open AI tests are commented out because OpenAI models are not being used in the prod.
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


def test_protecting_group_context_integration():
    """Test that protecting groups are properly detected and context is generated in LLM calls."""

    # Test molecule that contains protecting groups (benzyl ether pattern)
    molecule_with_pg = "COCc1ccccc1"  # Contains OBn protecting group pattern

    # Mock the LLM call to capture the prompt that gets sent
    with patch('src.utils.llm.litellm.completion') as mock_completion:
        # Set up mock response
        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[
            0].message.content = '{"data": [{"smiles": "test"}]}'
        mock_completion.return_value = mock_response

        # Call the LLM function
        status_code, response = call_LLM(molecule=molecule_with_pg,
                                         use_protecting_group_feature=True)

        # Verify the call was made
        assert status_code == 200
        assert mock_completion.called

        # Get the actual prompt that was sent to the LLM
        call_args = mock_completion.call_args
        messages = call_args[1]['messages']  # Get messages from kwargs
        user_message = None
        for msg in messages:
            if msg['role'] == 'user':
                user_message = msg['content']
                break

        # Verify protecting group context was included
        assert user_message is not None
        assert "PROTECTING GROUP CONTEXT" in user_message
        assert "Original SMILES: COCc1ccccc1" in user_message
        assert "Masked SMILES: %" in user_message
        assert "Symbol mapping:" in user_message
        assert "deprotection conditions:" in user_message


def test_no_protecting_group_context_when_disabled():
    """Test that protecting group context is NOT generated when feature is disabled."""

    molecule_with_pg = "COCc1ccccc1"  # Contains protecting group

    with patch('src.utils.llm.litellm.completion') as mock_completion:
        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[
            0].message.content = '{"data": [{"smiles": "test"}]}'
        mock_completion.return_value = mock_response

        # Call with protecting group feature disabled
        status_code, response = call_LLM(molecule=molecule_with_pg,
                                         use_protecting_group_feature=False)

        assert status_code == 200

        # Get the prompt
        call_args = mock_completion.call_args
        messages = call_args[1]['messages']
        user_message = None
        for msg in messages:
            if msg['role'] == 'user':
                user_message = msg['content']
                break

        # Verify NO protecting group context was included
        assert user_message is not None
        assert "PROTECTING GROUP CONTEXT" not in user_message


def test_no_protecting_group_context_for_non_pg_molecules():
    """Test that no protecting group context is generated for molecules without protecting groups."""

    # Simple molecule without protecting groups
    simple_molecule = "CC(=O)O"  # Acetic acid

    with patch('src.utils.llm.litellm.completion') as mock_completion:
        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[
            0].message.content = '{"data": [{"smiles": "test"}]}'
        mock_completion.return_value = mock_response

        status_code, response = call_LLM(molecule=simple_molecule,
                                         use_protecting_group_feature=True)

        assert status_code == 200

        # Get the prompt
        call_args = mock_completion.call_args
        messages = call_args[1]['messages']
        user_message = None
        for msg in messages:
            if msg['role'] == 'user':
                user_message = msg['content']
                break

        # Should not have protecting group context since no PGs detected
        assert user_message is not None
        assert "PROTECTING GROUP CONTEXT" not in user_message


@patch('src.utils.llm.get_config_info')
def test_protecting_group_context_uses_current_config(mock_get_config_info):
    """Test that protecting group context reflects the current configuration."""

    # Mock config info to return specific protecting group info
    mock_get_config_info.return_value = {
        'pg_count': 2,
        'pg_names': ['TestPG1', 'TestPG2']
    }

    molecule_with_pg = "CO"  # Simple methoxy that should be detected

    with patch('src.utils.llm.litellm.completion') as mock_completion:
        with patch('src.utils.llm.mask_protecting_groups_multisymbol'
                   ) as mock_mask:
            # Mock the masking function to show a protecting group was detected
            mock_mask.return_value = "$"  # Masked version

            mock_response = MagicMock()
            mock_response.choices = [MagicMock()]
            mock_response.choices[
                0].message.content = '{"data": [{"smiles": "test"}]}'
            mock_completion.return_value = mock_response

            status_code, response = call_LLM(molecule=molecule_with_pg,
                                             use_protecting_group_feature=True)

            assert status_code == 200

            # Verify that masking function was called
            mock_mask.assert_called_once_with(molecule_with_pg)

            # Get the prompt to verify context was generated
            call_args = mock_completion.call_args
            messages = call_args[1]['messages']
            user_message = None
            for msg in messages:
                if msg['role'] == 'user':
                    user_message = msg['content']
                    break

            # Should contain protecting group context
            assert user_message is not None
            assert "PROTECTING GROUP CONTEXT" in user_message


if __name__ == '__main__':
    pytest.main()
