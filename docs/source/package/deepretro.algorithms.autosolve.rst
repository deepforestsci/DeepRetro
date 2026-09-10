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

Route metadata
--------------

``solved`` is true when AZ returns a successful route for the target, or when
AZ successfully closes every terminal branch of an LLM-proposed route.
``az_solved`` is true only for a complete route whose reactions all came from
AZ. A completed mixed LLM/AZ route therefore has ``solved=True`` and
``az_solved=False``. A target AZ recognizes as available without reactions
has both flags true and an empty ``steps`` list.

Every reaction produced by autosolve carries ``solved_by`` with value ``"llm"``
or ``"az"``, both on the raw reaction node and on the parsed step. For example,
an LLM disconnection followed by an AZ reaction yields step sources
``["llm", "az"]``. The tag identifies the source of the disconnection; an LLM
step can still appear in an incomplete route. Terminal molecules are not
reaction steps. Agent proposals use the ``"llm"`` tag, including proposals
made after consulting AZ tools.

``az_summary.az_solved_all`` reports whether AZ closed every terminal branch,
so it can be true for mixed routes. The leaf counts remain available for
debugging. Older outputs used ``az_solved`` to mean that *any* leaf came from
AZ; consumers should use ``az_summary.leaves_az_generated > 0`` for that query.
Generic externally supplied trees with unknown provenance are not assigned
an inferred source.

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
