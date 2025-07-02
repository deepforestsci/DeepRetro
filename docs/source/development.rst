Development Guide
=================

Essential information for developing RecursiveLLM.

Setup
-----

.. code-block:: bash

   git clone <repository-url>
   cd recursiveLLM
   conda env create -f environment.yml
   conda activate dfs_si_challenge
   pip install -r tests/requirements_tests.txt

Project Structure
----------------

.. code-block:: text

   recursiveLLM/
   ├── src/                    # Main source code
   │   ├── api.py             # Flask API server (main entry point)
   │   ├── main.py            # Core retrosynthesis function
   │   ├── cache.py           # Caching functionality
   │   ├── prithvi.py         # Core retrosynthesis logic
   │   ├── rec_prithvi.py     # Recursive retrosynthesis
   │   ├── metadata.py        # Metadata extraction
   │   └── utils/             # Utility modules
   ├── viewer/                # Web interface
   ├── tests/                 # Test suite
   └── docs/                  # Documentation

Development Workflow
-------------------

1. **Create feature branch:**

   .. code-block:: bash

      git checkout -b feature/your-feature-name

2. **Make changes** following code standards
3. **Run tests:**

   .. code-block:: bash

      python -m pytest tests/

4. **Submit pull request**

Code Standards
--------------

* Follow PEP 8
* Use type hints
* Add docstrings (Google style)
* Keep functions focused

Example:

.. code-block:: python

   def process_molecule(smiles: str, model: str) -> Dict[str, Any]:
       """Process molecule with specified model.
       
       Args:
           smiles: SMILES string
           model: Model identifier
           
       Returns:
           Processing results
       """
       return result

Testing
-------

**Unit Tests:**

.. code-block:: python

   def test_parse_response():
       response = "test response"
       result = parse_response(response)
       assert result is not None

**Integration Tests:**

.. code-block:: python

   def test_retrosynthesis_api():
       response = client.post(
           '/api/retrosynthesis',
           headers={'X-API-KEY': 'test-key'},
           json={'smiles': 'CC'}
       )
       assert response.status_code == 200

**Mock External Services:**

.. code-block:: python

   @patch('src.utils.llm.call_LLM')
   def test_llm_integration(mock_llm):
       mock_llm.return_value = '{"result": "test"}'
       result = process_with_llm("test")
       assert result is not None

Adding Features
--------------

**New LLM Model:**

1. Add to `variables.py`:

   .. code-block:: python

      NEW_MODELS = ["new-model-name"]

2. Update `api.py` model selection
3. Add provider logic in `llm.py`

**New Validation:**

1. Create validation function
2. Add to main function with flag
3. Update API endpoint parameters

**New AiZynthFinder Model:**

1. Add to `AZ_MODEL_LIST` in `variables.py`
2. Update model validation in `api.py`

Error Handling
--------------

Use structured error handling:

.. code-block:: python

   try:
       result = process_data(input_data)
   except ValueError as e:
       return jsonify({"error": str(e)}), 400
   except Exception as e:
       logger.error(f"Unexpected error: {e}")
       return jsonify({"error": "Internal server error"}), 500

Logging
-------

Use structured logging:

.. code-block:: python

   import structlog
   
   log = structlog.get_logger()
   log.info("Processing molecule", smiles=smiles, model=model)

Documentation
-------------

* Update docstrings for new functions
* Update user guide for new features
* Update API documentation if needed
* Keep README current

For contribution guidelines, see :doc:`contributing`. 