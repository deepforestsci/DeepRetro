deepretro.utils.parse
=====================

Utilities for converting retrosynthesis route trees into the step and
dependency format consumed by the route viewer.

Overview
--------

``RetrosynthesisRouteParser`` is the primary API for new code. It keeps route
formatting state inside a small class and accepts injectable chemistry
callbacks, which makes the parser easier to test without loading heavyweight
chemistry dependencies.

The parser emits the viewer schema used by DeepRetro route visualizations:

- ``steps`` contains products, reactants, reagents, conditions, and reaction
  metrics for each parsed reaction step.
- ``dependencies`` maps each step id to the upstream step ids that produce its
  reactants.
- ``reactionmetrics`` contains ``scalabilityindex`` and
  ``closestliterature``.

The historical module-level functions remain available:

- ``parse_step`` parses a route tree into raw steps and dependencies.
- ``fix_dependencies`` rebuilds dependencies from product/reactant matches.
- ``format_output`` parses a route tree and returns viewer-ready output.

Example
-------

.. code-block:: python

   from deepretro.utils.parse import RetrosynthesisRouteParser

   parser = RetrosynthesisRouteParser(
       basic_molecules={"O"},
       chemical_formula_calculator=lambda smiles: smiles,
       mass_calculator=lambda smiles: 1.0,
       scalability_calculator=lambda reactant, product: "N/A",
   )
   output = parser.format_output(
       {
           "smiles": "CCOC(=O)N",
           "children": [
               {"children": [{"smiles": "CCO"}, {"smiles": "NC=O"}, {"smiles": "O"}]}
           ],
       }
   )

   assert output["steps"][0]["products"][0]["smiles"] == "CCOC(=O)N"
   assert [molecule["smiles"] for molecule in output["steps"][0]["reactants"]] == [
       "CCO",
       "NC=O",
   ]
   assert [molecule["smiles"] for molecule in output["steps"][0]["reagents"]] == ["O"]
   assert output["dependencies"] == {"1": []}

API Reference
-------------

.. automodule:: deepretro.utils.parse
   :members:
   :undoc-members:
