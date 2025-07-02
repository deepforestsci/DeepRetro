User Guide
==========

Complete guide for using RecursiveLLM effectively for retrosynthesis analysis.

Quick Start
-----------

Authentication
~~~~~~~~~~~~~~

All API requests require the ``X-API-KEY`` header:

.. code-block:: bash
   :linenos:
   :caption: Health check with authentication

   curl -H "X-API-KEY: your-api-key" \
        -H "Content-Type: application/json" \
        http://localhost:5000/api/health

Basic Retrosynthesis
~~~~~~~~~~~~~~~~~~~~

**Simple request:**

.. code-block:: bash
   :linenos:
   :caption: Basic retrosynthesis analysis

   curl -X POST http://localhost:5000/api/retrosynthesis \
     -H "X-API-KEY: your-key" \
     -H "Content-Type: application/json" \
     -d '{
       "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
       "model_type": "claude37",
       "advanced_prompt": true
     }'

**Python equivalent:**

.. code-block:: python
   :linenos:
   :caption: Python API request

   import requests

   response = requests.post(
       "http://localhost:5000/api/retrosynthesis",
       headers={"X-API-KEY": "your-key"},
       json={
           "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
           "model_type": "claude37",
           "advanced_prompt": True
       }
   )
   
   result = response.json()
   print(f"Found {len(result['data']['pathway'])} synthesis steps")

API Parameters
--------------

**Request Parameters:**

.. list-table::
   :widths: 20 15 15 15 35
   :header-rows: 1

   * - Parameter
     - Type
     - Required
     - Default
     - Description
   * - ``smiles``
     - string
     - ✓
     - 
     - Target molecule SMILES string
   * - ``model_type``
     - string
     - 
     - ``claude3``
     - LLM model: ``claude3``, ``claude37``, ``deepseek``, ``gpt4o``
   * - ``advanced_prompt``
     - boolean
     - 
     - ``false``
     - Enhanced prompting for better results
   * - ``model_version``
     - string
     - 
     - ``USPTO``
     - AiZynthFinder model version
   * - ``stability_flag``
     - boolean
     - 
     - ``false``
     - Enable molecular stability checks
   * - ``hallucination_check``
     - boolean
     - 
     - ``false``
     - Enable hallucination detection

Model Configuration
-------------------

LLM Models
~~~~~~~~~~

.. list-table::
   :widths: 25 15 20 40
   :header-rows: 1

   * - Model
     - Identifier
     - Best For
     - Internal Name
   * - **Claude 3 Opus**
     - ``claude3``
     - High quality results
     - ``claude-3-opus-20240229``
   * - **Claude 3.7 Sonnet**
     - ``claude37``
     - Hybrid reasoning
     - ``anthropic/claude-3-7-sonnet-20250219``
   * - **DeepSeek-R1**
     - ``deepseek``
     - Cost-effective analysis
     - ``fireworks_ai/accounts/fireworks/models/deepseek-r1``
   * - **GPT-4o**
     - ``gpt4o``
     - General purpose
     - ``gpt-4o``

.. note::
   **Claude 3.7 Sonnet** is the latest model with hybrid reasoning capabilities, 
   providing both fast responses and deep analytical thinking.

AiZynthFinder Models
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Model Version
     - Access
     - Description
   * - ``USPTO``
     - Free
     - Standard USPTO database (default)
   * - ``Pistachio_25``
     - Licensed
     - 25% Pistachio database coverage
   * - ``Pistachio_50``
     - Licensed  
     - 50% Pistachio database coverage
   * - ``Pistachio_100``
     - Licensed
     - 100% Pistachio database coverage
   * - ``Pistachio_100+``
     - Licensed
     - Enhanced Pistachio coverage with optimizations

Advanced Features
-----------------

Enhanced Prompting
~~~~~~~~~~~~~~~~~~

Enable advanced prompting for better results:

.. code-block:: python
   :linenos:
   :caption: Advanced prompting example

   response = requests.post(
       "http://localhost:5000/api/retrosynthesis",
       headers={"X-API-KEY": "your-key"},
       json={
           "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
           "model_type": "claude37",
           "advanced_prompt": True,           # Enhanced reasoning
           "stability_flag": True,            # Check stability
           "hallucination_check": True,       # Detect hallucinations
           "model_version": "Pistachio_50"    # Better database
       }
   )

Partial Rerun
~~~~~~~~~~~~~

Rerun analysis from a specific step:

.. code-block:: python
   :linenos:
   :caption: Partial rerun from specific step

   # First, get original analysis
   original = requests.post(
       "http://localhost:5000/api/retrosynthesis",
       headers={"X-API-KEY": "your-key"},
       json={"smiles": "target_molecule"}
   )
   
   # Then rerun from step 2 with different molecule
   rerun = requests.post(
       "http://localhost:5000/api/partial_rerun",
       headers={"X-API-KEY": "your-key"},
       json={
           "step_id": "step_2",
           "new_smiles": "CC(C)(C)OC(=O)Cl",
           "model_type": "claude37",
           "advanced_prompt": True
       }
   )

Complete Rerun
~~~~~~~~~~~~~~

Rerun entire analysis with updated parameters:

.. code-block:: python
   :linenos:
   :caption: Complete rerun with new settings

   response = requests.post(
       "http://localhost:5000/api/rerun_retrosynthesis",
       headers={"X-API-KEY": "your-key"},
       json={
           "model_type": "claude37",      # Switch to different model
           "advanced_prompt": True,       # Enable advanced features
           "stability_flag": True         # Add stability checks
       }
   )

Response Format
---------------

Success Response
~~~~~~~~~~~~~~~~

.. code-block:: json
   :linenos:
   :caption: Successful retrosynthesis response

   {
     "status": "success",
     "data": {
       "pathway": [
         {
           "step": 1,
           "step_id": "step_1", 
           "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
           "precursors": [
             {
               "smiles": "CC(C)(C)OC(=O)Cl",
               "confidence": 0.85,
               "reaction_type": "acylation",
               "availability": "commercial"
             },
             {
               "smiles": "N[C@@H](CC1=CC=CC=C1)C(=O)O",
               "confidence": 0.92,
               "reaction_type": "acylation", 
               "availability": "commercial"
             }
           ],
           "reaction_confidence": 0.88,
           "feasibility_score": 0.75
         }
       ],
       "metadata": {
         "model_used": "anthropic/claude-3-7-sonnet-20250219",
         "processing_time": 2.5,
         "total_steps": 1,
         "advanced_prompt": true,
         "stability_checked": true
       }
     }
   }

Error Response
~~~~~~~~~~~~~~

.. code-block:: json
   :linenos:
   :caption: Error response format

   {
     "status": "error",
     "error": {
       "code": "INVALID_SMILES",
       "message": "The provided SMILES string is invalid",
       "details": "Could not parse SMILES: 'CC(C)(C)OC(=O)N[C@@H]'"
     }
   }



Web Interface
-------------

The web interface provides an intuitive way to interact with RecursiveLLM:

Features
~~~~~~~~

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Feature
     - Description
   * - **Interactive Input**
     - Enter SMILES strings or paste from clipboard
   * - **Model Selection**
     - Choose from available LLM models
   * - **Pathway Visualization**
     - Interactive tree view of synthesis pathways
   * - **Confidence Indicators**
     - Visual confidence scores for each step
   * - **Step Editing**
     - Edit and rerun specific pathway steps
   * - **File Management**
     - Upload/download JSON pathway files
   * - **Export Options**
     - Export results as JSON, CSV, or images

Access
~~~~~~

Start the web interface:

.. code-block:: bash
   :linenos:
   :caption: Start web interface

   # Start backend
   python src/api.py
   
   # Start frontend (in new terminal)
   cd viewer
   python -m http.server 8000
   
   # Open browser
   # http://localhost:8000

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Invalid SMILES Error:**
- Validate SMILES strings using RDKit before submission
- Check for proper molecular notation format
- Ensure SMILES represents valid chemical structures

**API Key Issues:**
- Verify API_KEY environment variable is set
- Test with the `/api/health` endpoint to validate key
- Check for unauthorized (401) responses

**Model Availability:**
- Available models: claude3, claude37, deepseek
- AiZynthFinder models depend on local installation
- Check variables.py for current model list

Performance Optimization
------------------------

Caching Strategy
~~~~~~~~~~~~~~~~
- Results are automatically cached to improve performance
- Use `/api/clear_molecule_cache` to refresh cached data
- Cache persists across server restarts

Rate Limiting
~~~~~~~~~~~~~
- API implements rate limiting to prevent overload
- Add delays between requests for batch processing
- Handle 429 status codes with retry logic 