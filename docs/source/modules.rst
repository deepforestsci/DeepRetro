Python Package API
==================

API documentation for the ``deepretro`` Python package.

**Package structure:**

- **deepretro** - Main package; 
- **deepretro.featurizers** — Reaction-step featurizers (DeepChem-compatible).
- **deepretro.algorithms** — Recursive retrosynthesis (:func:`rec_run_DeepRetro`, :func:`single_run_DeepRetro`).
- **deepretro.utils** — LLM calls (:func:`call_LLM`, :func:`llm_pipeline`), AiZynthFinder (:func:`run_az`), and feature extraction.

.. toctree::
   :maxdepth: 4

   package/deepretro
