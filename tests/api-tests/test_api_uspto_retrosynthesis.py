import requests
import json
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import MOLECULE_1, BASE_URL, ENDPOINTS, X_API_KEY, DEEPSEEK_FIREWORKS_MODEL, USPTO_MODEL, CLAUDE_ADV_MODEL


def test_retrosynthesis_uspto_deepseek_m1p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
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


def test_retrosynthesis_uspto_deepseek_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
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


def test_retrosynthesis_uspto_claude_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "llm": CLAUDE_ADV_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert isinstance(response.json(), dict)
    assert ['dependencies', 'steps'] == list(response.json().keys())


def test_retrosynthesis_uspto_claude_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "llm": CLAUDE_ADV_MODEL,
        "model_version": USPTO_MODEL})

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
