import pytest
import json
import requests

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import MOLECULE_1, BASE_URL, ENDPOINTS, X_API_KEY, DEEPSEEK_FIREWORKS_MODEL, PISTACHIO_MODEL, USPTO_MODEL

def test_retrosynthesis_deepseek_pistachio_p1_success():
    """Test retrosynthesis endpoint with
    Model: DeepSeek
    Model Version: Pistachio
    Prompt: Advance.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "True",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_retrosynthesis_deepseek_pistachio_p0_success():
    """Test retrosynthesis endpoint with
    Model: DeepSeek
    Model Version: Pistachio
    Prompt: Basic.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "False",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_rerun_retro_deepseek_pistachio_p1_success():
    """Test rerun_retro endpoint with
    Model: DeepSeek
    Model Version: Pistachio
    Prompt: Advance.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['rerun_retro']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "True",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_rerun_retro_deepseek_pistachio_p0_success():
    """Test rerun_retro endpoint with
    Model: DeepSeek
    Model Version: Pistachio
    Prompt: Basic.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['rerun_retro']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "False",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_rerun_retro_deepseek_uspto_p1_success():
    """Test rerun_retro endpoint with
    Model: DeepSeek
    Model Version: USPTO
    Prompt: Advance.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['rerun_retro']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "True",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_rerun_retro_deepseek_uspto_p0_success():
    """Test rerun_retro endpoint with
    Model: DeepSeek
    Model Version: USPTO
    Prompt: Basic.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['rerun_retro']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "False",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_retrosynthesis_deepseek_uspto_p1_success():
    """Test retrosynthesis endpoint with
    Model: DeepSeek
    Model Version: USPTO
    Prompt: Advance.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "True",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_retrosynthesis_deepseek_uspto_p0_success():
    """Test retrosynthesis endpoint with
    Model: DeepSeek
    Model Version: USPTO
    Prompt: Basic.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "False",
        "hallucination_check": "False",
        "advanced_prompt": "False",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_stability_hallucination_claude_p0_success():
    """Test retrosynthesis endpoint with
    Model: DeepSeek
    Model Version: Pistachio
    Prompt: Advance.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "True",
        "hallucination_check": "True",
        "advanced_prompt": "False",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())

def test_stability_hallucination_claude_p1_success():
    """Test retrosynthesis endpoint with
    Model: DeepSeek
    Model Version: Pistachio
    Prompt: Advance.

    Asserts
    -------
    status_code : int
        200 for a successful request.
    response : dict
        The response should be a dictionary with keys ['dependencies', 'steps'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "stability_flag": "True",
        "hallucination_check": "True",
        "advanced_prompt": "True",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())

if __name__ == '__main__':
    pytest.main()
