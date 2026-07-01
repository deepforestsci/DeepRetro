deepretro.algorithms.autosolve
==============================

Recursive AiZynthFinder + LLM retrosynthesis solver. At each node AZ is tried
first; when it cannot solve a molecule the LLM (or, under
``solve_mode="single_step_agent"``, a tool-calling agent) proposes precursors,
which are recursively solved until AZ solves every leaf. ``run_llm`` is the
single owner of validity, stability, and hallucination filtering, so the agent
path can never bypass the safety checks.

.. code-block:: python

   from deepretro.algorithms.autosolve import AutoSolver

   # Non-agentic pipeline (default) with dependency-injected runners for testing:
   solver = AutoSolver(
       az_runner=lambda smiles, model: (False, []),
       llm_runner=lambda molecule, **kw: ([["OC(=O)c1ccccc1O"]], ["hydrolysis"], [0.9]),
       hallucination_mode="none",
   )
   output = solver.autosolve("CC(=O)Oc1ccccc1C(=O)O")

   # Tool-calling agent mode:
   agent_solver = AutoSolver(
       solve_mode="single_step_agent",
       tool_backend="sandbox",       # let the model write and run Python
   )

   # Top-K candidate routes for the batch runner:
   routes = agent_solver.solve_multiple("CC(=O)Oc1ccccc1C(=O)O", k=3)

Modes
-----

- ``solve_mode``: ``"pipeline"`` (default, non-agentic), ``"single_step_agent"``
  (per-molecule tool-calling agent), or ``"orchestrator"`` (reserved; raises
  ``NotImplementedError``).
- ``tool_backend``: ``"structured"`` (validity/stability/hallucination tools) or
  ``"sandbox"`` (adds a ``run_python`` code-execution tool). Ignored in
  ``pipeline`` mode.

API
---

.. automodule:: deepretro.algorithms.autosolve
   :members:
   :undoc-members:
