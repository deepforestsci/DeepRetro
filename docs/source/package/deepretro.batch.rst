deepretro.batch
===============

End-to-end batch retrosynthesis runner. Downloads a CSV from a public Google
Sheet, (optionally) trains the hallucination checker, reads target molecules
from a text file, runs each through
:class:`deepretro.algorithms.autosolve.AutoSolver`, and dumps the routes as
``<out>/<timestamp>/<molecule>/pathway_<i>.json``.

Output layout
-------------

::

   <out>/<timestamp>/<molecule-slug>/
       pathway_1.json       parsed, scored route
       pathway_2.json       ...
       llm_calls.jsonl      one JSON line per LLM call made for this target
       error.json           written instead of pathways when the molecule failed

``llm_calls.jsonl`` is described in :doc:`../logging` (section
"Per-molecule LLM call logs"); its records carry the Langfuse ``session_id`` that
groups the same calls server-side.

Command line
------------

.. code-block:: bash

   python scripts/run_batch.py \
       --sheet-url "https://docs.google.com/spreadsheets/d/<ID>/export?format=csv&gid=<GID>" \
       --molecules molecules.txt \
       --out batch_output \
       --solve-mode single_step_agent \
       --tool-backend sandbox \
       --agent-min-iterations 5 \
       --agent-max-iterations 15 \
       --agent-iteration-decay 0.75

``--sheet-url`` and ``--molecules`` are required. The three ``--agent-*``
flags size the per-node agent turn budget
(``clamp(carbons * decay**depth, min, max)``); the values shown are the
defaults. The training step is a
**template**: it runs only when the CSV carries ``product``/``reactants``/``label``
columns, otherwise the batch falls back to the heuristic hallucination checker.

API
---

.. automodule:: deepretro.batch
   :members:
   :undoc-members:
