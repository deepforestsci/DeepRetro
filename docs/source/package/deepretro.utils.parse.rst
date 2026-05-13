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

Input Tree Structure
--------------------

The retrosynthesis pipeline produces a recursive tree where each molecule
node may have ``children`` containing reaction wrappers, which in turn
contain precursor molecules:

.. code-block:: text

   root = {
       "smiles": "<product>",
       "children": [                    // reaction wrappers
           {
               "children": [            // precursor molecules
                   {"smiles": "<reactant_1>", "children": [...]},
                   {"smiles": "<reactant_2>"},
               ]
           }
       ]
   }

Algorithm
---------

The parser uses a depth-first traversal to convert this tree into a flat
list of reaction steps and a dependency map.

.. code-block:: text

   PARSE-NODE(node, S, D, parent_id)
   ──────────────────────────────────────────────────
   Input : node       — a route tree node
           S          — list of accumulated steps (mutated)
           D          — dependency map (mutated)
           parent_id  — step id of the calling parent, or NIL
   Output: S and D are updated in place
   ──────────────────────────────────────────────────
    1  step ← CREATE-STEP(node, |S| + 1)
    2  ATTACH-TO-PARENT(node, S, parent_id)
    3
    4  if step = NIL                          ▷ leaf node, no children
    5      if parent_id ≠ NIL
    6          D[parent_id] ← D[parent_id]    ▷ ensure key exists
    7      return
    8
    9  APPEND(S, step)
   10  if parent_id ≠ NIL
   11      APPEND(D[parent_id], step.id)
   12
   13  for each wrapper in node.children
   14      for each precursor in wrapper.children
   15          PARSE-NODE(precursor, S, D, step.id)

.. code-block:: text

   CREATE-STEP(node, step_id)
   ──────────────────────────────────────────────────
    1  if "children" ∉ node
    2      return NIL
    3  smiles ← node["smiles"]
    4  return {step: step_id, products: [smiles],
    5          reactants: [], reagents: [],
    6          reactionmetrics: [∅]}

.. code-block:: text

   ATTACH-TO-PARENT(node, S, parent_id)
   ──────────────────────────────────────────────────
    1  if parent_id = NIL or node.is_reaction
    2      return
    3  smiles ← node["smiles"]
    4  parent ← S[parent_id]
    5  if smiles ∈ basic_molecules
    6      APPEND(parent.reagents, smiles)
    7  else
    8      APPEND(parent.reactants, smiles)
    9  parent.scalability ← CALC-SCALABILITY(smiles, parent.product)

.. code-block:: text

   FORMAT-OUTPUT(root)
   ──────────────────────────────────────────────────
   Input : root — the root of a retrosynthesis route tree
   Output: {steps, dependencies}
   ──────────────────────────────────────────────────
    1  S ← [], D ← {}
    2  PARSE-NODE(root, S, D, NIL)
    3  D ← REBUILD-DEPENDENCIES(S)       ▷ overwrite tree-order deps
    4  return {steps: S, dependencies: D}

.. code-block:: text

   REBUILD-DEPENDENCIES(S)
   ──────────────────────────────────────────────────
    1  product_map ← {}
    2  for each step in S
    3      product_map[step.product.smiles] ← step.id
    4  D' ← {}
    5  for each step in S
    6      D'[step.id] ← []
    7      for each reactant in step.reactants
    8          if reactant.smiles ∈ product_map
    9              APPEND(D'[step.id], product_map[reactant.smiles])
   10  return D'

Example
-------

.. code-block:: python

   from deepretro.utils.parse import RetrosynthesisRouteParser

   parser = RetrosynthesisRouteParser(
       basic_molecules=set(),
       chemical_formula_calculator=lambda smiles: "N/A",
       mass_calculator=lambda smiles: 0.0,
       scalability_calculator=lambda reactant, product: "N/A",
   )
   output = parser.format_output(
       {
           "smiles": "CCO",
           "children": [{"children": [{"smiles": "CC"}, {"smiles": "O"}]}],
       }
   )

   assert output["steps"][0]["products"][0]["smiles"] == "CCO"
   assert output["dependencies"] == {"1": []}

API Reference
-------------

.. automodule:: deepretro.utils.parse
   :members:
   :undoc-members:
