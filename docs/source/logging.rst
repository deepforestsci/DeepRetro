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

API reference
-------------

.. currentmodule:: deepretro.logging

.. autofunction:: configure_logging
