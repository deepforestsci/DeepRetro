deepretro.utils.az
==================

AiZynthFinder integration helpers for template-based retrosynthesis.

What This Module Does
---------------------

- Runs AiZynthFinder search for a target SMILES.
- Returns route dictionaries with metadata/scores.
- Provides optional image outputs for generated routes.
- Short-circuits simple molecules to avoid unnecessary search overhead.

Configuration
-------------

The module resolves configuration in this order:

1. ``AZ_MODELS_PATH/<az_model>/config.yml`` (model-variant specific)
2. ``AZ_MODEL_CONFIG_PATH`` (fallback global config)

Both are interpreted relative to project root when loaded by this module.

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
