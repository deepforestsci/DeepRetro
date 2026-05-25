deepretro.algorithms.autosolve
==============================

``AutoSolver`` provides a package-level orchestration layer for single-molecule
retrosynthesis. It first attempts an AiZynthFinder route through
``deepretro.utils.az.run_az``. When no template route is available, it calls the
package LLM retrosynthesis pipeline, recursively solves returned precursors, and
formats the resulting tree with ``deepretro.utils.parse.format_output``.

Example
-------

.. code-block:: python

   from deepretro import AutoSolver

   solver = AutoSolver(
       llm="anthropic/claude-opus-4-6:adv",
       az_model="Pistachio_100+",
       hallucination_mode="heuristic",
   )
   result = solver.solve("CC(=O)Oc1ccccc1C(=O)O")
   print(result["solved"], len(result["steps"]))

API Reference
-------------

.. automodule:: deepretro.algorithms.autosolve
   :members:
   :undoc-members:
