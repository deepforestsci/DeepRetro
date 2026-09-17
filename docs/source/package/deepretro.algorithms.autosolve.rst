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

SMILES canonicalization
-----------------------

Every target SMILES is rewritten to its RDKit canonical form at the entry
point of ``solve``, ``single_step``, ``solve_multiple`` and ``autosolve``.
Model-proposed reactants are canonicalized in ``run_llm`` before the safety
filters. This means the recursion, the cycle-detection visited set, and the
emitted route tree all use one spelling per molecule. Unparseable strings are
left as-is and rejected by the validity filter. Batch input is also
canonicalized when the molecules file is read
(:func:`deepretro.batch.read_molecules`). Every system prompt tells the model
that its SMILES will be canonicalized. The ``autosolve`` output carries
``"smiles_canonicalized": true`` so downstream consumers know the route tree
uses canonical SMILES.

Modes
-----

- ``solve_mode``: ``"pipeline"`` (default, non-agentic), ``"single_step_agent"``
  (per-molecule tool-calling agent), or ``"orchestrator"`` (reserved; raises
  ``NotImplementedError``).
- ``tool_backend``: ``"structured"`` (validity/stability/hallucination tools) or
  ``"sandbox"`` (adds a ``run_python`` code-execution tool). Ignored in
  ``pipeline`` mode.

Agent iteration budget
----------------------

In ``single_step_agent`` mode each node gets its own turn budget instead of a
fixed ``max_iterations``. The budget is derived from the molecule's carbon
count and shrinks with recursion depth::

   budget = clamp(round(carbons * decay ** depth), min_iterations, max_iterations)

with defaults ``min_iterations=5``, ``max_iterations=15`` and ``decay=0.75``.
A 20-carbon target therefore gets 15 turns at the root, 11 at depth 2 and the
floor of 5 from depth 5 onward; a molecule RDKit cannot parse, or one without
carbon, gets the floor. Tune the three values through
``AutoSolver(agent_min_iterations=..., agent_max_iterations=...,
agent_iteration_decay=...)`` or the batch flags ``--agent-min-iterations``,
``--agent-max-iterations`` and ``--agent-iteration-decay``. The chosen budget
is logged per node as ``Agent iteration budget``. See
:func:`deepretro.agents.loop.iteration_budget`.

API
---

.. automodule:: deepretro.algorithms.autosolve
   :members:
   :undoc-members:
