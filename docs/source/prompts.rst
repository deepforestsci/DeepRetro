Prompt Strategy
===============

DeepRetro keeps its reusable prompt templates in
``deepretro/utils/variables.py``. The templates are organized by
provider-specific behavior, reasoning depth, and prompt modifiers so the
runtime can select the narrowest prompt that matches the requested model.

How prompt selection works
--------------------------

At runtime, ``src/utils/llm.py`` chooses a prompt pair based on the model name
and whether the model is invoked with the ``:adv`` suffix.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Runtime case
     - System prompt
     - User prompt
   * - Default providers
     - ``SYS_PROMPT``
     - ``USER_PROMPT``
   * - Default providers with ``:adv``
     - ``SYS_PROMPT_V4``
     - ``USER_PROMPT_V4``
   * - OpenAI-family models
     - ``SYS_PROMPT_OPENAI``
     - ``USER_PROMPT_OPENAI``
   * - DeepSeek-family models
     - ``SYS_PROMPT_DEEPSEEK``
     - ``USER_PROMPT_DEEPSEEK``
   * - DeepSeek-family models with ``:adv``
     - ``SYS_PROMPT_V4``
     - ``USER_PROMPT_DEEPSEEK_V4``

This split keeps the high-level retrosynthesis framing stable while allowing
provider-specific output constraints where the models need them.

Prompt families
---------------

The prompt constants fall into three groups:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Group
     - Purpose
   * - Baseline prompts
     - ``SYS_PROMPT`` and ``USER_PROMPT`` provide the default single-step
       retrosynthesis framing and JSON output contract.
   * - Advanced prompts
     - ``SYS_PROMPT_V4`` and ``USER_PROMPT_V4`` ask for a more explicit staged
       reasoning process and richer practical analysis.
   * - Provider-specific prompts
     - ``SYS_PROMPT_OPENAI``, ``USER_PROMPT_OPENAI``,
       ``SYS_PROMPT_DEEPSEEK``, ``USER_PROMPT_DEEPSEEK``, and
       ``USER_PROMPT_DEEPSEEK_V4`` adapt the output contract to known provider
       differences.

Prompt modifiers
----------------

The base prompt pair can be extended with targeted add-ons:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Constant
     - When it is added
   * - ``ADDON_PROMPT_7_MEMBER``
     - Appended when the molecule contains a seven-membered ring so the model
       pays extra attention to that motif.
   * - ``PROTECTING_GROUP_CONTEXT``
     - Appended when protecting groups are detected so the prompt explains the
       masking symbols and asks for real SMILES plus deprotection-aware
       reasoning.

Design intent
-------------

The prompt strategy follows a few simple rules:

1. Keep a stable baseline contract for ordinary retrosynthesis calls.
2. Isolate advanced reasoning prompts instead of mutating the baseline prompt in
   place.
3. Add provider-specific variants only where output formatting or model behavior
   diverges.
4. Layer structural context, such as protecting-group handling, as explicit
   modifiers instead of duplicating full prompt templates.

This structure makes prompt changes easier to review because each change tends
to fall into one of three categories: base behavior, provider tuning, or
context-specific add-ons.
