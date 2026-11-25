deepretro.utils.az
==================

AiZynthFinder integration helpers for template-based retrosynthesis.

Installation
------------

Install the package with ``deepretro[az]`` before using this module:

.. code-block:: bash

   uv pip install "deepretro[az]"

What This Module Does
---------------------

- Runs AiZynthFinder search for a target SMILES.
- Returns route dictionaries with metadata/scores.
- Provides optional image outputs for generated routes.
- Supports optional explicit caching through ``CacheManager``.
- Short-circuits simple molecules to avoid unnecessary search overhead.

Caching
-------

``run_az`` and ``run_az_with_img`` do not cache anything unless the caller
passes a cache instance explicitly:

.. code-block:: python

   from deepretro.utils.az import run_az
   from deepretro.utils.cache import CacheManager

   cache = CacheManager()
   solved, routes = run_az("C1CCCCC1", az_model="USPTO", cache=cache)

Configuration
-------------

The module resolves configuration in this order:

1. ``AZ_MODELS_PATH/<az_model>/config.yml`` (model-variant specific)
2. ``AZ_MODEL_CONFIG_PATH`` (fallback global config)

Path rules for both variables:

- An absolute path is used as given.
- A leading ``~`` is expanded to your home directory.
- A relative path is joined to ``PROJECT_ROOT``. This is the folder two
  levels above ``deepretro/utils/az.py``. In a source checkout, it is the
  repository root.
- When deepretro is installed as a package, ``PROJECT_ROOT`` is not the
  repository. Use absolute paths there.
- An unset variable is skipped.

The module reads both variables once, when it is imported. Set them before
you import ``deepretro.utils.az``.

If no model-specific config exists and ``AZ_MODEL_CONFIG_PATH`` is unset or
points to a missing file, the run raises ``FileNotFoundError``.

Basic Molecule Shortcut
-----------------------

``run_az`` and ``run_az_with_img`` bypass tree search if:

- The molecule is in ``BASIC_MOLECULES``, or
- ``is_basic_molecule(smiles)`` is ``True``

In these cases, the function returns a solved status and an in-stock molecule
route object immediately.


API
---

.. automodule:: deepretro.utils.az
   :members:
   :undoc-members:
