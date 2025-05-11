import requests
import json
import pytest

import rootutils
root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

from tests.variables_test import MOLECULE_1, BASE_URL, ENDPOINTS, X_API_KEY, DEEPSEEK_FIREWORKS_MODEL

def test_health_deepseek_success():
    """Test the health of the DeepSeek model with advanced_model=True and advanced_prompt=False.
    
    expected output:
        status_code: 200
        response: {'status': 'healthy'}
    """
  
    url = f"{BASE_URL}{ENDPOINTS['health']}"

    payload = json.dumps({
      "smiles": MOLECULE_1,
      "advanced_model": "True",
      "advanced_prompt": "False",
      "llm": DEEPSEEK_FIREWORKS_MODEL,
      "model_version": "USPTO"})

    headers = {
      'x-api-key': X_API_KEY,
      'Content-Type': 'application/json'
    }

    response = requests.request("GET", url, headers=headers, data=payload)

    assert response.status_code == 200
    assert response.json() == {'status': 'healthy'}


if __name__ == '__main__':
  pytest.main()
