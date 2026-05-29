deepretro.metadata
==================

LLM-assisted helpers for reaction metadata prediction.

The module exposes the ported metadata agents for reagent, condition, and
literature-reaction prediction. Functions return ``(status_code, payload)``
tuples for compatibility with the legacy metadata workflow and use package
imports, structured logging, and in-memory cache helpers.

``recommend_reaction_metadata`` is the high-level entry point for reaction
SMILES strings. It parses ``reactants>agents>products`` or
``reactants>>products`` input, recommends reagents, predicts reaction
conditions, and retrieves the nearest literature reaction.

Example
-------

.. code-block:: python

   from deepretro.metadata import reagent_agent

   status, reagents = reagent_agent(
       reactants=[{"smiles": "CCO"}],
       product=[{"smiles": "CC=O"}],
       model="openai/gpt-4o-mini",
   )

.. code-block:: python

   from deepretro.metadata import recommend_reaction_metadata

   status, metadata = recommend_reaction_metadata(
       "CCO.CC(=O)O>>CCOC(C)=O",
       model="openai/gpt-4o-mini",
   )

API Reference
-------------

.. automodule:: deepretro.metadata
   :members:
