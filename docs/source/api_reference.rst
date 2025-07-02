API Reference
=============

Essential API endpoints for retrosynthesis analysis.

Authentication
--------------

All requests require the ``X-API-KEY`` header.

Main Endpoints
--------------

- **POST /api/retrosynthesis**: Perform retrosynthesis analysis on a target molecule.
- **POST /api/partial_rerun**: Rerun retrosynthesis from a specific step with new parameters.
- **POST /api/rerun_retrosynthesis**: Rerun complete retrosynthesis with same or updated parameters.

Utility Endpoints
-----------------

- **GET /api/health**: Check API server health and status.
- **POST /api/clear_molecule_cache**: Clear cached results for a specific molecule.

Model Configuration
-------------------

**Available LLM Models:**

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Model
     - Identifier
     - Internal Name
   * - Claude 3 Opus
     - ``claude3``
     - ``claude-3-opus-20240229``
   * - Claude 3.7 Sonnet
     - ``claude37``
     - ``anthropic/claude-3-7-sonnet-20250219``
   * - DeepSeek-R1
     - ``deepseek``
     - ``fireworks_ai/accounts/fireworks/models/deepseek-r1``

**Available AiZynthFinder Models:**

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Model Version
     - Description
   * - ``USPTO``
     - Standard USPTO database (free, default)
   * - ``Pistachio_25``
     - 25% Pistachio database coverage
   * - ``Pistachio_50``
     - 50% Pistachio database coverage
   * - ``Pistachio_100``
     - 100% Pistachio database coverage
   * - ``Pistachio_100+``
     - Enhanced Pistachio coverage

Error Handling
--------------

**HTTP Status Codes:**

.. list-table::
   :widths: 15 25 60
   :header-rows: 1

   * - Code
     - Status
     - Description
   * - 200
     - OK
     - Request successful
   * - 400
     - Bad Request
     - Invalid parameters or SMILES
   * - 401
     - Unauthorized
     - Invalid or missing API key
   * - 404
     - Not Found
     - Endpoint not found
   * - 500
     - Internal Error
     - Server error

**Error Response Format:**

.. code-block:: json
   :linenos:
   :caption: Error response example

   {
     "status": "error",
     "error": {
       "code": "INVALID_SMILES",
       "message": "The provided SMILES string is invalid",
       "details": "Could not parse SMILES: 'invalid_smiles'"
     }
   }

**Common Error Codes:**

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Error Code
     - Description
   * - ``INVALID_SMILES``
     - SMILES string cannot be parsed
   * - ``UNAUTHORIZED``
     - API key is invalid or missing
   * - ``MODEL_NOT_FOUND``
     - Specified model is not available
   * - ``PROCESSING_FAILED``
     - Retrosynthesis analysis failed
   * - ``RATE_LIMIT_EXCEEDED``
     - Too many requests in time window

Rate Limiting
-------------

API requests are rate-limited to ensure fair usage:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Limit Type
     - Restriction
   * - Requests per minute
     - 60 requests/minute per API key
   * - Concurrent requests
     - 5 simultaneous requests per API key
   * - Daily requests
     - 10,000 requests/day per API key

When rate limits are exceeded, the API returns HTTP 429 with:

.. code-block:: json
   :linenos:

   {
     "status": "error",
     "error": {
       "code": "RATE_LIMIT_EXCEEDED",
       "message": "Rate limit exceeded. Please try again later.",
       "retry_after": 60
     }
   } 