deepretro.utils.visualize
=========================

Static pathway visualization for retrosynthesis outputs.

What This Module Does
---------------------

Renders the ``dict`` returned by
:meth:`deepretro.algorithms.autosolve.AutoSolver.solve` as a single
PIL image showing the full synthesis tree, mirroring the interactive
web viewer's ``processData`` + ``renderGraph`` behavior:

- Step 0 (the final target) is drawn as a large green circle on the left.
- Each subsequent step's reactants are drawn as blue circles, stacked
  vertically and connected to their parent via bezier edges.
- Layout traverses the full dependency tree left-to-right.

Output Image Layout
-------------------

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Visual element
     - Meaning
   * - Green circle
     - Step 0 target molecule (product of step 1).
   * - Blue circles
     - Reactants at each step, stacked vertically per column.
   * - Bezier curves
     - Dependency edges (parent step -> precursor step).
   * - Text above a circle
     - Step label, e.g. ``"Step 3"``.
   * - Text below a circle
     - ``<mass> g/mol`` and chemical formula.

Dependencies
------------

- **Pillow** (``PIL``) for canvas construction and text.
- **RDKit** for SMILES parsing and 2D structure rendering (with an
  automatic fall-back from the Cairo drawer to ``Draw.MolToImage``
  when Cairo is unavailable).
- **NumPy** for the alpha-channel transparency pass that blends
  structures over the coloured node circles.

Molecule metadata (formula, mass) is read from the solver output when
present (``product_metadata`` / ``reactant_metadata``) and otherwise
recomputed on the fly via RDKit.

Font Handling
-------------

``get_font`` tries the following TrueType files in order and falls back
to Pillow's built-in bitmap font if none are available:

1. ``arial.ttf``
2. ``Arial.ttf``
3. ``DejaVuSans.ttf``

This keeps the module usable in minimal container images without
bundled fonts.

Example
-------

.. code-block:: python

   from deepretro.algorithms.autosolve import AutoSolver
   from deepretro.utils.visualize import visualize_pathway

   solver = AutoSolver()
   result = solver.solve("CC(=O)Oc1ccccc1C(=O)O")
   img = visualize_pathway(result)
   img.save("aspirin_pathway.png")

Empty results (``{"steps": [], "dependencies": {}}``) return a small
placeholder image so calling code never has to branch on emptiness.


API
---

.. automodule:: deepretro.utils.visualize
   :members:
   :undoc-members:
