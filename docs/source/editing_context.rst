Editing Context
===============

This page complements :doc:`knowledge_graph` with the minimum file context you
should load before editing a hotspot.

Runtime Flow Hotspots
---------------------

.. list-table::
   :widths: 24 28 28 20
   :header-rows: 1

   * - If you edit
     - Read first
     - Why
     - Watch for
   * - ``src/api.py``
     - ``config/advanced_settings.json``, ``src/main.py``, ``viewer/config.js``
     - Owns request validation, model flags, rerun endpoints, and the contract
       the browser calls
     - Changing payloads can break reruns and frontend defaults
   * - ``src/main.py``
     - ``src/prithvi.py``, ``src/utils/custom_logging.py``
     - Thin entry point into the runtime orchestration path
     - Logging setup changes affect all request execution
   * - ``src/prithvi.py``
     - ``src/rec_prithvi.py``, ``src/utils/parse.py``, ``src/metadata.py``
     - Bridges recursion, route formatting, and metadata enrichment
     - Output shape changes surface immediately in the viewer
   * - ``src/rec_prithvi.py``
     - ``src/utils/az.py``, ``src/utils/llm.py``, ``src/utils/job_context.py``
     - Controls template-first fallback logic and recursive precursor expansion
     - Small branching changes can alter the whole synthesis search behavior
   * - ``src/utils/llm.py``
     - ``src/variables.py``, ``src/utils/stability_checks.py``,
       ``src/utils/hallucination_checks.py``
     - Owns prompt selection, provider routing, response parsing, and
       validation filters
     - Model identifiers and prompt constants are tightly coupled
   * - ``src/utils/parse.py``
     - ``src/utils/utils_molecule.py``, ``viewer/index.html``
     - Converts route trees into viewer-consumable dependency graphs
     - Dependency rebuilding is sensitive to duplicate product SMILES
   * - ``src/metadata.py``
     - ``src/variables.py``, ``src/cache.py``
     - Adds reagents, conditions, and literature metadata after route creation
     - Extra LLM calls and cache behavior affect latency noticeably

Package Hotspots
----------------

.. list-table::
   :widths: 24 28 28 20
   :header-rows: 1

   * - If you edit
     - Read first
     - Why
     - Watch for
   * - ``deepretro/__init__.py``
     - ``deepretro/featurizers/reactionstep.py``, ``deepretro/logging.py``,
       ``deepretro/models/hallucination_classifier.py``
     - Defines the intended public import surface
     - New exports change how consumers discover the package API
   * - ``deepretro/featurizers/reactionstep.py``
     - ``deepretro/utils/domain_features.py``, ``deepretro/tests/test_reactionstep.py``
     - Core ML representation for product/reactant reaction steps
     - Feature dimension changes ripple into loaders and models
   * - ``deepretro/data/loader.py``
     - ``deepretro/featurizers/reactionstep.py``, ``deepretro/tests/test_loader.py``
     - Streams reaction CSV rows into DeepChem datasets
     - Column assumptions and dropped-row behavior matter to training code
   * - ``deepretro/models/hallucination_classifier.py``
     - ``deepretro/featurizers/reactionstep.py``, ``deepretro/utils/metrics.py``
     - Wraps the featurizer and threshold logic into a trainable classifier
     - Model persistence and threshold behavior are easy places to regress
   * - ``deepretro/utils/az.py``
     - ``src/cache.py``, ``src/variables.py``, ``deepretro/tests/test_az.py``
     - Package-facing AiZynthFinder helper with legacy runtime coupling
     - Refactors here can silently break package/runtime boundary assumptions
   * - ``deepretro/logging.py``
     - ``deepretro/tests/test_logging.py``
     - Canonical structlog configuration for package consumers
     - Contextvar behavior is covered by focused tests and should stay stable

Frontend Hotspots
-----------------

.. list-table::
   :widths: 24 28 28 20
   :header-rows: 1

   * - If you edit
     - Read first
     - Why
     - Watch for
   * - ``viewer/index.html``
     - ``viewer/config.js``, ``src/api.py``, ``config/advanced_settings.json``
     - Owns API calls, model selection UI, and partial rerun interactions
     - Request payload drift breaks the main user path immediately
   * - ``viewer/config.js``
     - ``src/api.py``, ``config/advanced_settings.json``
     - Maps frontend endpoint names and default request options
     - Model key mismatches cause confusing UI/runtime divergence

Non-Obvious Dependencies
------------------------

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Dependency
     - Why it matters
   * - ``deepretro/utils/az.py -> src/cache.py``
     - The package AiZynthFinder helper still uses the runtime cache decorator,
       so package-only refactors can reach into the runtime tree.
   * - ``deepretro/utils/az.py -> src/variables.py``
     - Basic-molecule heuristics come from runtime constants instead of a pure
       package-local definition.
   * - ``viewer/index.html -> config/advanced_settings.json``
     - The viewer reads model defaults that must stay aligned with backend
       validation logic.
   * - ``src/utils/parse.py -> viewer/index.html``
     - The parser shapes the route/dependency graph that the viewer expects to
       render and rerun.
   * - ``src/api.py -> partial rerun state``
     - Full rerun and partial rerun logic share request/response assumptions
       that are easy to break when only testing one endpoint.

Suggested Review Order
----------------------

1. Load the graph on :doc:`knowledge_graph`.
2. Open the hotspot row for the file you want to change.
3. Read the listed context files before editing.
4. Run the targeted tests for that cluster instead of only broad smoke checks.
