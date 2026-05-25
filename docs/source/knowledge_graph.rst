Knowledge Graph
===============

This page is the fastest way to orient yourself in the DeepRetro repository.
It is organized around the paths that control runtime behavior, reusable
package code, and the browser UI.

.. raw:: html

   <div style="margin: 1rem 0 1.5rem 0;">
     <iframe
       src="_static/knowledge_graph.html"
       title="DeepRetro knowledge graph"
       style="width: 100%; height: 860px; border: 1px solid #d7dce5; border-radius: 12px; background: #ffffff;"
     ></iframe>
   </div>

If the embedded graph does not load, open the standalone version at
``docs/source/_static/knowledge_graph.html`` in the built docs output.

How To Read The Graph
---------------------

- **Runtime nodes** are the production request path under ``src/``.
- **Package nodes** are the reusable ``deepretro/`` modules intended for
  notebooks, training, and library-style imports.
- **Frontend nodes** are the browser entry points under ``viewer/``.
- **Config nodes** are files that control runtime model choices and defaults.
- **Cross-links** highlight hidden or non-obvious dependencies, especially when
  a package module still imports from the runtime tree.

Start Here
----------

If you need to make a change quickly, use these starting points:

- API or request validation: ``src/api.py``
- Recursive retrosynthesis flow: ``src/prithvi.py`` and ``src/rec_prithvi.py``
- LLM prompts, parsing, or validation gates: ``src/utils/llm.py``
- Viewer JSON shaping: ``src/utils/parse.py``
- UI behavior and reruns: ``viewer/index.html`` and ``viewer/config.js``
- Reusable ML featurization: ``deepretro/featurizers/reactionstep.py``
- DeepChem dataset ingestion: ``deepretro/data/loader.py``
- Package hallucination model: ``deepretro/models/hallucination_classifier.py``

Critical Repo Facts
-------------------

- ``src/`` is the active Flask runtime. It is what the web app uses today.
- ``deepretro/`` is a separate reusable package. It is not the live API path.
- The repository has split packaging metadata: the root ``pyproject.toml``
  packages from ``src``, while ``deepretro/pyproject.toml`` packages
  ``deepretro*`` from the repository root.
- ``deepretro/utils/az.py`` still depends on ``src.cache`` and
  ``src.variables``. That is the most important cross-boundary coupling to know
  before refactoring.

Next Step
---------

For file-by-file editing context, read :doc:`editing_context`.
