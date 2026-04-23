deepretro.algorithms.hallucination\_checker
============================================

Heuristic hallucination checker for retrosynthetic reaction steps.
Compares reactant and product molecules to detect structural
inconsistencies (atom-count mismatches, ring-size changes, substituent
swaps, etc.) and produces a 0–100 hallucination score.

.. code-block:: python

   from deepretro.algorithms.hallucination_checker import (
       hallucination_compare_molecules,
       calculate_hallucination_score,
   )

   issues = hallucination_compare_molecules("c1ccccc1", "c1ccccc1OC")
   result = calculate_hallucination_score("c1ccccc1", "c1ccccc1OC")
   print(result["score"], result["severity"])

API
---

.. automodule:: deepretro.algorithms.hallucination_checker
   :members:
   :undoc-members:
