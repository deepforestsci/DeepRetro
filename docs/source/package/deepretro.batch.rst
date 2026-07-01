deepretro.batch
===============

End-to-end batch retrosynthesis runner. Downloads a CSV from a public Google
Sheet, (optionally) trains the hallucination checker, reads target molecules
from a text file, runs each through
:class:`deepretro.algorithms.autosolve.AutoSolver`, and dumps the routes as
``<out>/<timestamp>/<molecule>/pathway_<i>.json``.

Command line
------------

.. code-block:: bash

   python scripts/run_batch.py \
       --sheet-url "https://docs.google.com/spreadsheets/d/<ID>/export?format=csv&gid=<GID>" \
       --molecules molecules.txt \
       --out batch_output \
       --solve-mode single_step_agent \
       --tool-backend sandbox

``--sheet-url`` and ``--molecules`` are required. The training step is a
**template**: it runs only when the CSV carries ``product``/``reactants``/``label``
columns, otherwise the batch falls back to the heuristic hallucination checker.

API
---

.. automodule:: deepretro.batch
   :members:
   :undoc-members:
