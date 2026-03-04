Stability Checker
=================

Heuristic stability analysis for molecules proposed by retrosynthesis models.

.. code-block:: python

   from deepretro.algorithms import check_molecule_stability

   result = check_molecule_stability("c1ccccc1")
   print(result["assessment"])       # "Likely stable"
   print(result["stability_score"])  # 100

.. automodule:: deepretro.algorithms.stability_checker
   :members:
   :undoc-members:
