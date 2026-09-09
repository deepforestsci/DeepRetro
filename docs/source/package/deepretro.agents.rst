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

Protecting-group abstraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both tool backends expose ``handle_protection`` and ``handle_deprotection``.
The first hides existing candidate Boc, Cbz, TBS, OBn, OEt, and OMe motifs
behind mapped dummy atoms so the LLM can focus on core disconnections. The
second supports two explicit modes:

* ``mode="restore"`` (default): restore hidden atoms using private, per-agent
  restoration data. Requires ``mask_id``; this does not change chemical identity.
* ``mode="propose"``: propose chemical deprotection products from full protected
  SMILES. Omit ``mask_id``; an optional ``groups`` list restricts motifs.

.. code-block:: python

   registry = build_tool_registry()
   masked = registry.execute("handle_protection", {
       "smiles": "COc1ccc(C=O)cc1", "groups": ["OMe"]
   })
   # The agent reasons about the exposed core while retaining each [*:n].
   restored = registry.execute("handle_deprotection", {
       "smiles": masked["masked_smiles"], "mask_id": masked["mask_id"]
   })
   assert restored["smiles"] == masked["original_smiles"]

   proposed = registry.execute("handle_deprotection", {
       "smiles": "CCNC(=O)OC(C)(C)C", "mode": "propose", "groups": ["Boc"]
   })
   assert proposed["candidates"][0]["product_smiles"] == "CCN"

Proposal mode returns one candidate per recognized site, removing only that
motif and adding a hydrogen to the retained O/N atom. Other groups, scaffold
stereochemistry, atom maps, and salt components remain intact. Site indices
refer to canonical ``protected_smiles``; ``site_atom_map`` preserves the input
map number (zero when unmapped). The output explicitly reports the **forward**
direction, full substrate/product structures, and ``status="structural_candidate"``.
For a retrosynthetic deprotection step, supply a proposed protected precursor
and compare the resulting product with the target; do not use a deprotected
product as a retro precursor of its protected input.

These proposals are structural hypotheses, with ``conditions=None``. Conditions,
selectivity (including whether other similar sites also react), and substrate
compatibility require assessment. ``reaction_smiles`` omit reagents and
byproducts and are not atom-balanced. Empty candidates mean no supported site
was recognized, not chemical impossibility. Proposal mode rejects masked input,
does not require an earlier masking call, and does not change restoration state.

Pass each masked precursor separately to restoration using the same
``mask_id``; each may contain a subset of the original markers. IDs are valid
only in the registry that created them. Unknown or duplicate markers and
altered attachment types return tool errors. Full structures retain atom maps,
stereochemistry, and disconnected components. The agent is instructed to
restore before validation; final answers containing dummy atoms trigger a
repair turn within the existing iteration budget.

Matches identify structural motifs, not their intended synthetic role. The
optional ``groups`` list restricts matching for masking and proposal mode
(an empty list disables it); restoration rejects that argument.
Unsupported, bridging, multiply attached motifs, and hidden fragments with
specified atom stereochemistry remain explicit. Chemical proposals also exclude
stereogenic O/N attachment atoms. The
tools do not recommend installation reagents, validate reaction compatibility,
or guarantee that the model preserves every hidden group across a proposed
step; the agent must assess those questions using the full structures.

This adapts the masking-and-prompt idea from the history of
``protect_star_update`` (before ``cc4b536``). It uses graph operations instead
of that version's lossy text substitutions. The branch's later protection-site
recommendation subsystem is separate from these tools.

.. automodule:: deepretro.agents.protection
   :members:

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
