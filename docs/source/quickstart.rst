Quickstart Guide
================

Get up and running with DeepRetro in 5 minutes.

Installation
------------

**Option 1: Docker Installation (Recommended)**

.. code-block:: bash
   :linenos:
   :caption: Docker setup

   git clone <repository-url>
   cd recursiveLLM
   cp env.example .env
   # Edit .env with your API keys
   docker-compose up -d

**Option 2: Local Development Installation**

**1. Clone and setup environment:**

.. code-block:: bash
   :linenos:
   :caption: Environment setup

   git clone <repository-url>
   cd DeepRetro
   conda env create -f environment.yml
   conda activate deepretro

**2. Configure API keys:**

Create `.env` file:

.. code-block:: bash
   :linenos:
   :caption: Environment configuration

   # Required: Backend API key
   API_KEY=your-backend-key
   
   # LLM API keys (choose based on your model preference)
   ANTHROPIC_API_KEY=your-anthropic-key      # For Claude models
   FIREWORKS_API_KEY=your-fireworks-key      # For DeepSeek models

**3. Download models:**

.. code-block:: bash
   :linenos:
   :caption: Model download

   mkdir -p aizynthfinder/models
   python -c "
   from aizynthfinder.utils.download_public_data import download_public_data
   download_public_data('aizynthfinder/models/')
   "

Quick Start
-----------

**1. Start the backend:**

.. code-block:: bash
   :linenos:
   :caption: Start API server

   python src/api.py

**2. Start the web interface:**

.. code-block:: bash
   :linenos:
   :caption: Start web interface

   cd viewer
   python -m http.server 8000

.. important::

   Before using the web interface, **edit** ``viewer/config.js`` **to set your backend API endpoint** (e.g., http://localhost:5000). This ensures the frontend communicates with your running backend server.

**3. Test the API:**

.. code-block:: bash
   :linenos:
   :caption: Test API endpoint

   curl -H "X-API-KEY: your-api-key" \
        -H "Content-Type: application/json" \
        http://localhost:5000/api/health

Available Models
================

**LLM Models:**

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Model
     - Identifier
   * - **Claude 3 Opus**
     - ``claude3``
   * - **Claude 3.7 Sonnet**
     - ``claude37``
   * - **Claude 4 Sonnet**
     - ``claude4``
   * - **DeepSeek-R1**
     - ``deepseek``

**AiZynthFinder Models:**

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Model
     - Access
     - Description
   * - ``USPTO``
     - Free
     - Standard USPTO database (default, downloaded automatically in Docker)
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
     - Enhanced coverage with optimizations

Basic Usage
-----------

API Request
~~~~~~~~~~~

**Simple analysis:**

.. code-block:: bash
   :linenos:
   :caption: Basic API request

   curl -X POST http://localhost:5000/api/retrosynthesis \
     -H "X-API-KEY: your-key" \
     -H "Content-Type: application/json" \
     -d '{
       "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
       "model_type": "claude37"
     }'

**Advanced analysis:**

.. code-block:: bash
   :linenos:
   :caption: Advanced features enabled

   curl -X POST http://localhost:5000/api/retrosynthesis \
     -H "X-API-KEY: your-key" \
     -H "Content-Type: application/json" \
     -d '{
       "smiles": "CC(C)(C)OC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)O",
       "model_type": "claude37",
       "advanced_prompt": true,
       "stability_flag": true,
       "hallucination_check": true,
       "model_version": "USPTO"
     }'

Python Usage
~~~~~~~~~~~~

Make API requests using the `requests` library. See :doc:`api_reference` for complete endpoint documentation.

Web Interface
~~~~~~~~~~~~~

Open `http://localhost:8000` in your browser.

**Features:**

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Feature
     - Description
   * - **SMILES Input**
     - Enter SMILES strings or paste from clipboard
   * - **Model Selection**
     - Choose from Claude 3, Claude 3.7, Claude 4, DeepSeek
   * - **Interactive Visualization**
     - Tree view of synthesis pathways with confidence scores
   * - **Step Editing**
     - Edit and rerun specific pathway steps
   * - **File Management**
     - Upload/download JSON pathway files
   * - **Export Options**
     - Export as JSON, CSV, or images

Response Format
---------------

**Success Response:**

.. code-block:: json
   :linenos:
   :caption: Successful API response

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
         "advanced_prompt": true
       }
     }
   }

**Error Response:**

.. code-block:: json
   :linenos:
   :caption: Error response format

   {
     "status": "error",
     "error": {
       "code": "INVALID_SMILES",
       "message": "The provided SMILES string is invalid",
       "details": "Could not parse SMILES: 'invalid_string'"
     }
   }

Common Parameters
-----------------

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
     - —
     - Target molecule SMILES string
   * - ``model_type``
     - string
     - 
     - ``claude3``
     - LLM model: ``claude3``, ``claude37``, ``deepseek``, ``claude4opus``
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

Next Steps
----------

.. tip::
   **Ready to explore more?**
   
   * :doc:`user_guide` - Complete API documentation
   * :doc:`api_reference` - Detailed API reference
   * :doc:`development` - Development setup and contribution guide

**Common Use Cases:**

1. **Drug Discovery** - Analyze pharmaceutical intermediates
2. **Chemical Synthesis** - Plan multi-step organic syntheses  
3. **Process Development** - Optimize synthetic routes
4. **Research** - Explore novel synthetic pathways
5. **Education** - Learn retrosynthetic analysis

**Getting Help:**

* Check the :doc:`user_guide` for troubleshooting
* Open an issue for bugs or feature requests 