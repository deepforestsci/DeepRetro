Installation Guide
=================

Quick setup for DeepRetro.

Prerequisites
-------------

**Option 1: Docker (Recommended)**
* Docker and Docker Compose
* Git

**Option 2: Local Development**
* Python 3.9+
* Conda/Miniconda
* Git

System Requirements
------------------

* **OS**: Windows 10+, macOS 10.14+, Ubuntu 18.04+
* **Memory**: 16GB recommended
* **Storage**: 5GB free space

Installation Steps
-----------------

**Option 1: Docker Installation (Recommended)**

1. **Clone repository:**

   .. code-block:: bash

      git clone <repository-url>
      cd recursiveLLM

2. **Set up environment:**

   .. code-block:: bash

      cp env.example .env
      # Edit .env with your API keys

3. **Start the service:**

   .. code-block:: bash

      docker-compose up -d

4. **Verify installation:**

   .. code-block:: bash

      curl -H "X-API-KEY: your-api-key" http://localhost:5000/api/health

**Option 2: Local Development Installation**

1. **Clone repository:**

   .. code-block:: bash

      git clone <repository-url>
      cd DeepRetro

2. **Create environment:**

   .. code-block:: bash

      conda env create -f environment.yml
      conda activate deepretro

3. **Download models:**

   .. code-block:: bash

      mkdir -p aizynthfinder/models
      python -c "from aizynthfinder.utils.download_public_data import download_public_data; download_public_data('aizynthfinder/models/')"

4. **Configure API keys:**

   Create `.env` file:

   .. code-block:: bash

      API_KEY=your-backend-key
      ANTHROPIC_API_KEY=your-anthropic-key  # For Claude
      FIREWORKS_API_KEY=your-fireworks-key  # For DeepSeek

5. **Verify installation:**

   .. code-block:: bash

      python -c "import src.api; print('Installation successful!')"


Configuration
-------------

Model Configuration
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - LLM Model
     - Identifier
   * - Claude 3 Opus
     - ``claude3``
   * - Claude 3.7 Sonnet
     - ``claude37``
   * - Claude 4 Sonnet
     - ``claude4``
   * - DeepSeek-R1
     - ``deepseek``

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - AiZynthFinder Model
     - Description
   * - USPTO
     - Standard database (free, default, downloaded automatically in Docker)
   * - Pistachio_25
     - 25% Pistachio database (licensed)
   * - Pistachio_50
     - 50% Pistachio database (licensed)
   * - Pistachio_100
     - 100% Pistachio database (licensed)
   * - Pistachio_100+
     - Enhanced Pistachio coverage (licensed)

Getting Help
------------

* Check existing GitHub issues
* Create new issue with error details
* Review :doc:`user_guide` for troubleshooting 