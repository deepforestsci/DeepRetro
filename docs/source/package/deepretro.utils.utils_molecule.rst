deepretro.utils.utils_molecule
==============================

Molecule utilities for SMILES validation, substructure matching, molecular properties, and ring detection.

Overview
--------

The ``utils_molecule`` module provides chemistry-focused helpers used throughout the retrosynthesis pipeline:

- **SMILES validation** — Check validity and compare molecules
- **Substructure matching** — Query whether one molecule is a substructure of another
- **Molecular properties** — Weight, formula, fingerprints
- **Validity checks** — Filter LLM-proposed pathways for chemical validity and
  reject target-matching fragments
- **Ring detection** — Detect 7- and 8-member rings in molecules

Function Overview
-----------------

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Function
     - Purpose
   * - ``is_valid_smiles``
     - Check if a SMILES string parses to a valid molecule.
   * - ``substructure_matching``
     - Return 1 if query is a substructure of target, 0 otherwise.
   * - ``are_molecules_same``
     - Compare two SMILES (canonical form or fingerprint).
   * - ``validity_check``
     - Filter LLM pathways: keep valid precursors, drop same-as-target or substructures.
   * - ``calc_mol_wt``
     - Molecular weight from SMILES (returns 0.0 on invalid input).
   * - ``calc_chemical_formula``
     - Molecular formula from SMILES (returns "N/A" on invalid input).
   * - ``compute_fingerprint``
     - Morgan fingerprint as a bit vector list.
   * - ``detect_seven_member_rings``
     - True if molecule contains a 7-member ring.
   * - ``detect_eight_member_rings``
     - True if molecule contains an 8-member ring.


Usage
-----

.. code-block:: python

   from deepretro.utils.utils_molecule import (
       is_valid_smiles,
       substructure_matching,
       validity_check,
       calc_mol_wt,
       calc_chemical_formula,
       detect_seven_member_rings,
   )

   # Validate SMILES
   assert is_valid_smiles("CCO") is True
   assert is_valid_smiles("invalid!!!") is False

   # Substructure check (benzene in ethylbenzene)
   assert substructure_matching("CCc1ccccc1", "c1ccccc1") == 1

   # Filter LLM pathways
   pathways, explanations, confidence = validity_check(
       molecule="c1ccccc1",
       res_molecules=[["CC(=O)O", "c1ccccc1O"]],
       res_explanations=["ester hydrolysis"],
       res_confidence=[0.8],
   )

   # Molecular properties
   assert calc_mol_wt("CCO") > 0
   assert calc_chemical_formula("C") == "CH4"

   # Ring detection
   assert detect_seven_member_rings("C1CCCCCC1") is True
   assert detect_seven_member_rings("C1CCCCC1") is False

API
---

.. automodule:: deepretro.utils.utils_molecule
   :members:
   :undoc-members:
