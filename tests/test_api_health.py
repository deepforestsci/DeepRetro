import requests
import json
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import MOLECULE_1, BASE_URL, ENDPOINTS, X_API_KEY, DEEPSEEK_MODEL, CLAUDE_MODEL

# Health

# Deepseek Health
def test_health_deepseek_m1p1_success():
  
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
      "smiles": MOLECULE_1,
      "advanced_model": "True",
      "advanced_prompt": "True",
      "model_version": DEEPSEEK_MODEL})

    headers = {
      'x-api-key': X_API_KEY,
      'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_deepseek_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": DEEPSEEK_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


# Claude Health
def test_health_claude_m1p1_success():
  
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
      "smiles": MOLECULE_1,
      "advanced_model": "True",
      "advanced_prompt": "True",
      "model_version": CLAUDE_MODEL})

    headers = {
      'x-api-key': X_API_KEY,
      'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m0p1_success():
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "True",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m1p0_success():
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "True",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


def test_health_claude_m0p0_success():
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
        "smiles": MOLECULE_1,
        "advanced_model": "False",
        "advanced_prompt": "False",
        "model_version": CLAUDE_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200


if __name__ == '__main__':
  pytest.main()
