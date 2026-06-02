deepretro.algorithms.autosolve
==============================

``AutoSolver`` provides the main retrosynthesis pipeline for the deepretro
package. It combines AiZynthFinder template matching with LLM-based fallback,
route parsing, and metadata enrichment in a single class.

Pipeline
--------

.. code-block:: text

   autosolve(smiles)
     ├── solve(smiles)         → recursive AZ + LLM retrosynthesis → route tree
     ├── parse(route_tree)     → format into steps + dependencies
     └── add_metadata(parsed)  → enrich with reagents, conditions, literature

Individual methods can be called separately for finer control. Use
``single_step()`` for non-recursive one-pass retrosynthesis.

Example
-------

.. code-block:: python

   from deepretro import AutoSolver

   solver = AutoSolver(
       llm="anthropic/claude-opus-4-6:adv",
       az_model="Pistachio_100+",
       hallucination_mode="heuristic",
   )

   # Full pipeline
   result = solver.autosolve("CC(=O)Oc1ccccc1C(=O)O")
   print(result["solved"], len(result["steps"]))

   # Individual steps
   route_tree, solved = solver.solve("CC(=O)Oc1ccccc1C(=O)O")
   parsed = solver.parse(route_tree, solved=solved)
   enriched = solver.add_metadata(parsed)

   # Single-step (no recursion)
   route_tree, solved = solver.single_step("CC(=O)Oc1ccccc1C(=O)O")

API Reference
-------------

.. automodule:: deepretro.algorithms.autosolve
   :members:
   :undoc-members:
