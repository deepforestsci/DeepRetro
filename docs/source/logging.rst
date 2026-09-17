Logging Guide
=============

``deepretro`` uses `structlog <https://www.structlog.org/>`_ for structured,
leveled logging.  Every module creates a logger at the top of the file and emits
log events with key-value context that is machine-parseable and
human-readable.

Quick start
-----------

No configuration is required.  By default, structlog prints colored log lines to
stdout at all levels.

To control the log level or output format, call
:func:`deepretro.logging.configure_logging` once at your application entry
point:

.. code-block:: python

   from deepretro.logging import configure_logging

   # Show only warnings and above
   configure_logging(level="WARNING")

   # Silence all deepretro logging
   configure_logging(level="CRITICAL")

   # Verbose JSON output for production / log aggregation
   configure_logging(level="DEBUG", json_output=True)

If you already call ``structlog.configure()`` in your application (e.g. in
``src/main.py``), that configuration applies automatically — you do not need
``configure_logging``.

Minimal application example
---------------------------

Use ``configure_logging()`` once near your entry point, then create module
loggers with ``structlog.get_logger()``:

.. code-block:: python

   import structlog

   from deepretro.logging import configure_logging

   configure_logging(level="INFO")

   logger = structlog.get_logger(__name__)
   logger.info("Starting retrosynthesis", molecule="CCO")

For production log pipelines, switch the same helper to JSON output:

.. code-block:: python

   configure_logging(level="DEBUG", json_output=True)

How to log in a new module
--------------------------

Add two lines at the top of your module:

.. code-block:: python

   import structlog

   logger = structlog.get_logger()

Then use the standard log levels:

.. code-block:: python

   logger.debug("Low-level detail", key="value")
   logger.info("Normal operation", molecule=smiles)
   logger.warning("Recoverable issue", attempt=attempt)
   logger.error("Failure", smiles=smiles, error=str(exc))

Prefer **structured key-value pairs** over f-string interpolation so that
downstream consumers (ELK, Datadog, etc.) can parse fields automatically.

Log level conventions
---------------------

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Level
     - When to use
   * - ``DEBUG``
     - Detailed diagnostics: prompt text, raw response bodies, intermediate
       computation results.  Hidden by default in production.
   * - ``INFO``
     - Normal milestones: "Starting retrosynthesis", "Pipeline succeeded",
       "Validity check complete".
   * - ``WARNING``
     - Recoverable problems: retries, fallbacks, rejected pathways, missing
       optional config.
   * - ``ERROR``
     - Failures that prevent the current operation from completing: invalid
       SMILES, JSON parse errors, all API attempts exhausted.

Context propagation with ``bound_contextvars``
----------------------------------------------

Use ``structlog.contextvars.bound_contextvars`` to attach context (e.g. a
job ID) that will appear in **every** log line emitted within the block — even
from deeply nested functions:

.. code-block:: python

   from structlog.contextvars import bound_contextvars

   def run_retrosynthesis(molecule: str, job_id: str):
       with bound_contextvars(job_id=job_id, molecule=molecule):
           logger.info("Starting retrosynthesis")
           # All downstream log calls automatically include job_id + molecule
           result = call_llm(molecule)
           return result

This replaces the previous ``contextvars.ContextVar`` / ``job_context`` pattern.
There is no need to pass loggers through function arguments.

Per-molecule LLM call logs
--------------------------

Structured logs record *what the solver did*; the per-molecule LLM call log
records *what was sent to and returned by the model*. Both are written for every
target molecule the batch runner solves.

Where the file lives
~~~~~~~~~~~~~~~~~~~~

:func:`deepretro.batch.run_batch` opens a trace around each molecule, so the log
lands next to that molecule's routes::

   <out>/<timestamp>/<molecule-slug>/
       pathway_1.json
       pathway_2.json
       llm_calls.jsonl      <-- one JSON object per line, one line per LLM call
       error.json           (only when the molecule failed)

Outside the batch runner, pass ``llm_log_dir`` to
:class:`deepretro.algorithms.autosolve.AutoSolver` and
:meth:`~deepretro.algorithms.autosolve.AutoSolver.autosolve` writes to
``<llm_log_dir>/<molecule-slug>/llm_calls.jsonl``:

.. code-block:: python

   from deepretro.algorithms.autosolve import AutoSolver

   solver = AutoSolver(llm_log_dir="llm_logs")
   solver.autosolve("CC(=O)Oc1ccccc1C(=O)O")

Nothing is written when neither is configured, and a write failure is logged
once as a warning and then ignored — logging never breaks a run.

Record fields
~~~~~~~~~~~~~

Each line is one JSON object. ``kind`` says which of the two record shapes it
is; the first nine fields are shared by both:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Meaning
   * - ``timestamp``
     - ISO 8601 UTC time the record was written.
   * - ``kind``
     - ``llm_call`` (one model call) or ``tool_results`` (the tools the agent
       ran after one model turn).
   * - ``session_id``
     - Langfuse session id shared by every record of this run.
   * - ``trace_id``
     - Langfuse trace id shared by every generation and tool event of this run.
   * - ``target``
     - Root molecule the run is solving.
   * - ``node_molecule``
     - Molecule of the tree node that made the call.
   * - ``depth``
     - Retrosynthesis recursion depth of that node (``0`` at the root).
   * - ``stage``
     - ``retrosynthesis`` (LLM pipeline), ``retrosynthesis_agent`` (agent loop),
       or ``metadata`` (reagent/conditions/literature enrichment).
   * - ``iteration``
     - Agent loop turn or provider retry attempt, 1-based.

``llm_call`` records add:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Meaning
   * - ``model``
     - LiteLLM model identifier.
   * - ``messages``
     - Conversation sent to the provider. Opaque provider replay data
       (Anthropic ``thinking_blocks`` signatures, ``provider_specific_fields``)
       is stripped and tool calls are stored as plain JSON objects.
   * - ``response``
     - Assistant content string, or the serialized message when the turn only
       requested tools.
   * - ``tool_calls``
     - Tool calls the assistant requested, if any, as JSON objects.
   * - ``usage``
     - ``prompt_tokens`` / ``completion_tokens`` / ``total_tokens`` when the
       provider reports them.
   * - ``latency_ms``
     - Wall-clock duration of the provider call.
   * - ``error``
     - Error text when the call failed, otherwise ``null``.

``tool_results`` records add:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Meaning
   * - ``tool_results``
     - One object per executed tool: ``tool_call_id``, ``name``, ``arguments``
       (parsed) and ``output`` (the tool's return value). Written after the
       tools run, so the outputs of the agent's last allowed turn are kept
       even when no further model call follows.

Reading a log back is ordinary JSON Lines:

.. code-block:: python

   import json
   from pathlib import Path

   path = Path("batch_output/2026-07-01_00-00-00/CCO_1b0ef7e2/llm_calls.jsonl")
   records = [json.loads(line) for line in path.read_text().splitlines()]
   calls = [r for r in records if r["kind"] == "llm_call"]
   total_tokens = sum((r["usage"] or {}).get("total_tokens", 0) for r in calls)
   tool_outputs = [
       t for r in records if r["kind"] == "tool_results" for t in r["tool_results"]
   ]

How Langfuse traces are grouped
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every record of one run carries the same ``trace_id`` and ``session_id``.
LiteLLM receives both in the completion ``metadata``, so all of the run's
generations land in one Langfuse trace, each tagged with its ``stage``,
``depth`` and ``node_molecule`` as generation metadata. When
``LANGFUSE_PUBLIC_KEY`` and ``LANGFUSE_SECRET_KEY`` are set, every tool the
agent runs is also sent as a Langfuse *event* in that trace
(``name="tool:<tool name>"``, ``input`` = the tool arguments, ``output`` = the
tool result), so tool outputs render online beside the model turns that
requested them. Without the keys, only the local file is written.


LiteLLM's Langfuse callback sees one completion at a time, so DeepRetro supplies
the grouping itself. :func:`deepretro.utils.llm_trace.langfuse_metadata` adds the
keys LiteLLM's Langfuse integration understands to every completion's
``metadata``:

* ``session_id`` — the same value for every call of one target molecule, so the
  whole retrosynthesis appears as one Langfuse session. Generated as
  ``autosolve_<molecule-slug>_<YYYYmmdd_HHMMSS>_<pid>``.
* ``trace_name`` — always ``autosolve``.
* ``generation_name`` — the ``stage`` of the call.
* ``trace_metadata`` — ``molecule``, ``node_molecule``, and ``depth``.
* ``tags`` — ``deepretro``, ``autosolve``, and the stage.

Langfuse itself is enabled by LiteLLM through ``LANGFUSE_SECRET_KEY``,
``LANGFUSE_PUBLIC_KEY``, and ``LANGFUSE_HOST``. The local ``llm_calls.jsonl`` is
independent of it and needs no account or network access.

Traces nest by reuse, never by replacement: the batch runner opens the outer
trace per molecule and ``AutoSolver.autosolve`` reuses it, so one run never
splits across two sessions or two log files.

API reference
-------------

.. currentmodule:: deepretro.logging

.. autofunction:: configure_logging

.. currentmodule:: deepretro.utils.llm_trace

.. autoclass:: MoleculeTrace
   :members:

.. autofunction:: molecule_trace

.. autofunction:: current_trace

.. autofunction:: current_depth

.. autofunction:: current_node_molecule

.. autofunction:: node_depth

.. autofunction:: langfuse_metadata

.. autofunction:: record_llm_call
