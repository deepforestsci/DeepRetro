import requests
import json
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import BASE_URL, ENDPOINTS, X_API_KEY, DEEPSEEK_FIREWORKS_MODEL, USPTO_MODEL, CLAUDE_ADV_MODEL, PISTACHIO_MODEL

def test_retrosynthesis_fail():
    """Tests retrosynthesis endpoint with empty input("")

    Expected output:
        status_code: 500
        error message: status code and error message
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "advanced_model": "True",
        "advanced_prompt": "True",
        "llm": CLAUDE_ADV_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 400
    assert response.json() == {"error": "SMILES string is required. Please include a 'smiles' field"}


def test_rerun_retro_fail():
    """Tests rerun_retrosynthesis endpoint with empty input("")

    Expected output:
        status_code: 500
        error message: status code and error message
    """
    url = f"{BASE_URL}{ENDPOINTS['rerun_retro']}"

    payload = json.dumps({
        "advanced_model": "True",
        "advanced_prompt": "True",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 400
    assert response.json() == {"error": "Molecule string is required, Please include a 'smiles' field"}

if __name__ == '__main__':
    pytest.main()
