deepretro.utils.llm
===================

LiteLLM-backed utilities for single-step retrosynthesis prompts, completions,
response parsing, and pathway filtering.

Overview
--------

The module exposes a compact public API:

- ``obtain_prompt()`` selects the system/user prompt pair for Claude,
  DeepSeek, OpenAI, or generic model families.
- ``call_LLM()`` calls a real LiteLLM completion endpoint and returns
  ``(status_code, response_text)``.
- ``parse_response()`` extracts the JSON payload returned by Claude-style,
  DeepSeek-style, or OpenAI-style model output.
- ``validate_json_response()`` converts the parsed payload into pathways,
  explanations, and confidence scores.
- ``llm_pipeline()`` runs the call, parse, validation, and pathway filtering
  flow.

OpenAI Support
--------------

OpenAI models can be passed either with a LiteLLM provider prefix, such as
``openai/gpt-4o-mini``, or as OpenAI model names listed in
``deepretro.utils.variables.OPENAI_MODELS``. OpenAI chat models use
``max_completion_tokens`` and deterministic seed support when available.
OpenAI reasoning models, such as ``gpt-5`` or ``o``-series models, receive
``reasoning_effort`` when reasoning controls are enabled.

Example
-------

.. code-block:: python

   from deepretro.utils.llm import call_LLM, parse_response

   status, response_text = call_LLM(
       molecule="CCO",
       model="openai/gpt-4o-mini",
       max_output_tokens=2048,
   )

   if status == 200:
       parse_status, thinking_steps, json_content = parse_response(
           response_text,
           "openai/gpt-4o-mini",
       )

Testing
-------

The LLM endpoint test uses a real OpenAI call through LiteLLM and does not mock
the model. Set ``OPENAI_API_KEY`` in the environment or ``.env`` before running:

.. code-block:: bash

   uv run --project deepretro --extra dev python -m pytest deepretro/tests/test_llm.py -q

API
---

.. automodule:: deepretro.utils.llm
   :members:

.. automodule:: deepretro.utils.llm_helpers
   :members:
   :no-index:
