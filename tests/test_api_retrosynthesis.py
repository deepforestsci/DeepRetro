import requests
import json
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import MOLECULE_1, BASE_URL, ENDPOINTS, X_API_KEY, OPENAI_MODEL, DEEPSEEK_MODEL, CLAUDE_MODEL, EMPTY_RESPONSE, PISTACHIO_MODEL


def test_health_openai_m1p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_openai_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_openai_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_openai_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m1p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m1p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


# -------------------------------------------------------- #


def test_health_openai_m1p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_openai_m1p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_openai_m0p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_openai_m0p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": OPENAI_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_deepseek_m1p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_deepseek_m1p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_deepseek_m0p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_deepseek_m0p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_claude_m1p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_claude_m1p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_claude_m0p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_claude_m0p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


# Pistachio
def test_health_pistachio_m1p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_pistachio_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_pistachio_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_pistachio_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 200



def test_health_pistachio_m1p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_pistachio_m1p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_pistachio_m0p1_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500


def test_health_pistachio_m0p0_fail():
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "smiles": EMPTY_RESPONSE,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 500

if __name__ == '__main__':
    pytest.main()
