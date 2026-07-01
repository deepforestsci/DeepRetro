deepretro.agents
================

The agentic tool-calling layer used by
:class:`deepretro.algorithms.autosolve.AutoSolver` when
``solve_mode="single_step_agent"``. The model may call tools to validate SMILES,
check stability, check for hallucinations, or run Python in a sandbox, then emits
the same tagged-JSON payload the non-agentic pipeline produces.

.. code-block:: python

   from deepretro.agents.tools import build_tool_registry
   from deepretro.agents.sandbox import SubprocessSandbox
   from deepretro.agents.loop import agentic_single_step

   registry = build_tool_registry(tool_backend="sandbox", sandbox=SubprocessSandbox())
   registry.execute("validate_smiles", {"smiles": "C(O)C"})  # -> canonical 'CCO'

Sandbox security
----------------

``SubprocessSandbox`` runs model-written code in a subprocess that is
**secret-scrubbed** (no provider API keys inherited), resource-limited
(``RLIMIT_AS`` / ``RLIMIT_CPU`` / ``RLIMIT_NPROC`` / ``RLIMIT_FSIZE``),
wall-clock bounded, and — on Linux with ``unshare`` — network-isolated. It is
adequate for *semi-trusted* model-generated code; untrusted data is never
interpolated into executed code. For stronger isolation, implement the
``Sandbox`` protocol with a container backend (e.g. Podman ``--network=none``).

Tools API
---------

.. automodule:: deepretro.agents.tools
   :members:
   :undoc-members:

Sandbox API
-----------

.. automodule:: deepretro.agents.sandbox
   :members:
   :undoc-members:

Agent loop API
--------------

.. automodule:: deepretro.agents.loop
   :members:
   :undoc-members:
