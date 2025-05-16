import requests
import json
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import BASE_URL, ENDPOINTS, X_API_KEY, \
    DEEPSEEK_FIREWORKS_MODEL, USPTO_MODEL, CLAUDE_MODEL, PISTACHIO_MODEL


def test_retrosynthesis_fail():
    """Tests retrosynthesis endpoint with empty input("")

    Asserts
    -------
    status_code : int
            400 for a failed request.
    error message : dict
            The response should be a dictionary with keys ['error'].
    """
    url = f"{BASE_URL}{ENDPOINTS['retrosynthesis']}"

    payload = json.dumps({
        "advanced_prompt": "True",
        "stability_flag": "False",
        "hallucination_check": "False",
        "llm": CLAUDE_MODEL,
        "model_version": PISTACHIO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 400
    assert response.json() == {
        "error": "SMILES string is required. Please include a 'smiles' field"}


def test_rerun_retro_fail():
    """Tests rerun_retrosynthesis endpoint with empty input("")

    Asserts
    -------
        status_code: 400
        error message: status code and error message
    """
    url = f"{BASE_URL}{ENDPOINTS['rerun_retro']}"

    payload = json.dumps({
        "advanced_prompt": "True",
        "stability_flag": "False",
        "hallucination_check": "False",
        "llm": DEEPSEEK_FIREWORKS_MODEL,
        "model_version": USPTO_MODEL})

    headers = {
        'x-api-key': X_API_KEY,
        'Content-Type': 'application/json'
    }

    response = requests.request("POST", url, headers=headers, data=payload)

    assert response.status_code == 400
    assert response.json() == {
        "error":
        "Molecule string is required, Please include a 'smiles' field"}


if __name__ == '__main__':
    pytest.main()
