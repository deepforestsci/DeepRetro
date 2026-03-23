deepretro.utils
===============

Utility layer for domain features and template-based retrosynthesis integration.

This module group supports two major workflows:

- Feature engineering utilities for ML-ready reaction-step vectors.
- AiZynthFinder orchestration helpers for template route generation.

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


AiZynthFinder Integration Notes
-------------------------------

``run_az`` and ``run_az_with_img`` use environment-configured model paths:

- ``AZ_MODELS_PATH`` (preferred model-variant path)
- ``AZ_MODEL_CONFIG_PATH`` (fallback config path)

Behavior highlights:

- Auto-bypass for trivial/basic molecules via ``BASIC_MOLECULES`` and
  ``is_basic_molecule``.
- Caching via ``src.cache.cache_results`` decorator.
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

API Reference
-------------

.. automodule:: deepretro.utils
   :members:

deepretro.utils.domain\features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: deepretro.utils.domain_features
   :members:

