deepretro.utils
===============

Utility layer for domain features and template-based retrosynthesis integration.

This module group supports three major workflows:

- Feature engineering utilities for ML-ready reaction-step vectors.
- AiZynthFinder orchestration helpers for template route generation.
- Molecule utilities for SMILES validation, substructure matching, and pathway filtering (``deepretro.utils.utils_molecule``).

Utility Overview
----------------

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Utility
     - Purpose
   * - ``run_az``
     - Execute AiZynthFinder and return solved flag + route dictionaries.
   * - ``run_az_with_img``
     - Same as ``run_az`` plus route images when available.
   * - ``is_basic_molecule``
     - Shortcut heuristic for trivial molecules to bypass heavy search.
   * - ``CacheManager`` / ``make_cache_key``
     - Explicit in-memory cache primitives for expensive library operations.
   * - ``call_LLM`` / ``llm_pipeline``
     - LiteLLM-backed retrosynthesis calls, response parsing, and pathway filtering.


AiZynthFinder Integration Notes
-------------------------------

``run_az`` and ``run_az_with_img`` use environment-configured model paths:

- ``AZ_MODELS_PATH`` (preferred model-variant path)
- ``AZ_MODEL_CONFIG_PATH`` (fallback config path)

Behavior highlights:

- Auto-bypass for trivial/basic molecules via ``BASIC_MOLECULES`` and
  ``is_basic_molecule``.
- Caching via ``src.cache.cache_results`` decorator.
- Explicit process-local caching helpers live in ``deepretro.utils.cache`` when
  callers need in-memory caching without shared global state.
- Returns route dictionaries with metadata and scores from AiZynthFinder.

Example: run template search
----------------------------

.. code-block:: python

   from deepretro.utils.az import run_az

   solved, routes = run_az("C1CCCCC1", az_model="USPTO")
   print(solved, len(routes))


Submodules
----------

.. toctree::
   :maxdepth: 1

   deepretro.utils.az
   deepretro.utils.metrics
   deepretro.utils.cache
   deepretro.utils.llm
   deepretro.utils.utils_molecule

API Reference
-------------

.. automodule:: deepretro.utils
   :members:

deepretro.utils.domain_features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: deepretro.utils.domain_features
   :members:
