deepretro.utils.llm
===================

LiteLLM-backed single-step retrosynthesis utilities.

The LLM code is split into three layers:

- ``deepretro.utils.llm`` is the public workflow facade. It exposes prompt
  selection, model calls, response parsing, JSON validation, pathway filtering,
  and the ``llm_pipeline()`` orchestration function.
- ``deepretro.utils.llm_interface`` owns provider-specific behavior. The
  ``LLMInterface`` base class defines prompt construction, completion parameter
  construction, API calling, and response parsing. Concrete implementations
  handle Claude-style, OpenAI-style, DeepSeek-style, and generic responses.
- ``deepretro.utils.llm_helpers`` contains model normalization and small parsing
  helpers shared by the interface layer.

Provider Interfaces
-------------------

Use ``create_llm_interface()`` when code needs direct access to the
provider-specific interface:

.. code-block:: python

   from deepretro.utils.llm_interface import LLMRequest, create_llm_interface

   interface = create_llm_interface("openai/gpt-4o-mini")
   request = LLMRequest(
       molecule="CCO",
       model="openai/gpt-4o-mini",
       max_output_tokens=2048,
       enable_thinking=False,
   )

   messages = interface.build_messages(request)
   params = interface.build_completion_params(request, messages)

The factory returns one of these implementations:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Interface
     - Provider family
     - Parser behavior
   * - ``AnthropicLLM``
     - Anthropic / Claude
     - Requires ``<cot>`` content with at least one ``<thinking>`` entry and a
       JSON payload. JSON is accepted in ``<json>`` tags or as fenced/raw JSON
       after ``</cot>`` because Claude can return either shape.
   * - ``OpenAILLM``
     - OpenAI
     - Extracts tagged, fenced, or raw JSON and does not return thinking steps.
   * - ``DeepSeekLLM``
     - DeepSeek
     - Extracts optional ``<think>`` content plus tagged, fenced, or raw JSON.
   * - ``GenericLLM``
     - Fallback
     - Uses the Claude-style parser for compatible providers.

Public Workflow
---------------

Most callers should use the facade functions in ``deepretro.utils.llm``:

.. code-block:: python

   from deepretro.utils.llm import call_LLM, parse_response, validate_split_json

   status, response_text = call_LLM(
       molecule="CCO",
       model="openai/gpt-4o-mini",
       max_output_tokens=2048,
       enable_thinking=False,
   )

   if status == 200:
       parse_status, thinking_steps, json_content = parse_response(
           response_text,
           "openai/gpt-4o-mini",
       )
       if parse_status == 200:
           validate_split_json(json_content)

``llm_pipeline()`` combines the same steps into the end-to-end flow:

.. code-block:: python

   from deepretro.utils.llm import llm_pipeline

   pathways, explanations, confidence = llm_pipeline(
       molecule="CCO",
       model="openai/gpt-4o-mini",
       stability_check=False,
       hallucination_check=False,
       max_output_tokens=2048,
       enable_thinking=False,
   )

Model Selection
---------------

Model identifiers are normalized before calling LiteLLM:

- OpenAI models can be passed with a LiteLLM prefix, for example
  ``openai/gpt-4o-mini``, or as names recognized by
  ``deepretro.utils.variables.OPENAI_MODELS``.
- OpenAI chat models use ``max_completion_tokens`` and receive a deterministic
  ``seed``.
- OpenAI reasoning models, such as ``gpt-5`` and ``o``-series models, use
  ``reasoning_effort`` when reasoning controls are enabled. Their output-token
  budget is raised to a provider-safe minimum.
- Anthropic Claude 4 Opus/Sonnet models also receive ``reasoning_effort`` when
  reasoning controls are enabled. Their output-token budget is raised to the
  same provider-safe minimum, and their temperature is set to ``1`` for
  reasoning calls.
- DeepSeek aliases such as ``fireworks/deepseek-v3p2`` are normalized to the
  preferred Fireworks-hosted DeepSeek R1 model.
- A ``:adv`` suffix, such as ``openai/gpt-4o-mini:adv``, selects the advanced
  prompt mode unless an explicit ``prompt_mode`` argument is provided.

Testing
-------

Most tests in this module do not require live LLM credentials. They exercise
provider selection, completion-parameter construction, parser behavior, JSON
validation, and pipeline orchestration using local test doubles.

The focused test file also includes two slow Anthropic integration tests. When
``ANTHROPIC_API_KEY`` is configured, they call the live server with the real
retrosynthesis prompt for aspirin in both standard and advanced prompt modes.
Those tests assert that the response contains ``<cot>`` / ``<thinking>`` tags,
that JSON can be extracted from either tagged or fenced output, and that
``validate_split_json()`` can convert the payload into aligned pathways,
explanations, and confidence scores.

Run the focused test file with:

.. code-block:: bash

   uv run --project deepretro pytest deepretro/tests/test_llm.py -q

Run only the live Anthropic prompt checks with:

.. code-block:: bash

   ANTHROPIC_API_KEY=... uv run --project deepretro pytest \
       deepretro/tests/test_llm.py::test_live_anthropic_retrosynthesis_prompt_returns_parseable_tagged_json -q

API
---

.. automodule:: deepretro.utils.llm
   :members:

.. automodule:: deepretro.utils.llm_interface
   :members:
   :no-index:

.. automodule:: deepretro.utils.llm_helpers
   :members:
   :no-index:
