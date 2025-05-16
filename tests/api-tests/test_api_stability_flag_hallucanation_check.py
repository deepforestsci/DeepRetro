import pytest
import json
import requests

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import MOLECULE_1, BASE_URL, ENDPOINTS, X_API_KEY, USPTO_MODEL, CLAUDE_MODEL, PISTACHIO_MODEL


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
        "llm": CLAUDE_MODEL,
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
        "llm": CLAUDE_MODEL,
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