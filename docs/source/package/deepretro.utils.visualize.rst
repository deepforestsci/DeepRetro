deepretro.utils.visualize
=========================

Static pathway visualization for retrosynthesis outputs.

What This Module Does
---------------------

Renders the viewer-ready pathway format used by DeepRetro as a single
PIL image showing the full synthesis tree:

- Step 0 (the final target) is drawn as a large green circle on the left.
- Each subsequent step usually renders reactants as blue circles, stacked
  vertically and connected to their parent via bezier edges. Terminal
  product-only steps fall back to products so they remain visible.
- Layout traverses the full dependency tree left-to-right.

Scope
-----

This module renders retrosynthesis outputs as static PIL images.
SVG export is intentionally out of scope for this PR.

Input Shape
-----------

The renderer expects the same viewer-ready route shape used elsewhere in
the package:

- ``steps``: a flat list of reaction steps
- ``dependencies``: a mapping from step ids to precursor-producing step ids

At minimum, step ``1`` should exist and should contain a product molecule.

Dependencies
------------

- **Pillow** (``PIL``) for canvas construction and text.
- **RDKit** for SMILES parsing and 2D structure rendering.
- **NumPy** for the alpha-channel transparency pass that blends
  structures over the colored node circles.

Pillow is treated as an optional visualization dependency.
Install it with ``pip install "deepretro[viz]"``.

Example
-------

.. code-block:: python

   from deepretro.utils.visualize import visualize_pathway

   result = {
       "steps": [
           {
               "step": "1",
               "reactants": [{"smiles": "CCO"}],
               "products": [{"smiles": "CC=O"}],
           }
       ],
       "dependencies": {"1": []},
   }

   img = visualize_pathway(result)
   img.save("pathway.png")

Empty results (``{"steps": [], "dependencies": {}}``) return a small
placeholder image so calling code never has to branch on emptiness.

API
---

.. automodule:: deepretro.utils.visualize
   :members:
   :undoc-members:
