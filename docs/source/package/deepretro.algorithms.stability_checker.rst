Stability Checker
=================

Heuristic stability checker for molecules in retrosynthetic pathways.

When an LLM proposes reactant molecules, some of them may be chemically
unstable — strained small rings, anti-aromatic systems, reactive
intermediates like carbocations or carbenes, etc.  This module catches
those problems by inspecting the molecular graph with RDKit descriptors
and SMARTS pattern matching.  No trained model is required.

.. code-block:: python

   from deepretro.algorithms import check_molecule_stability

   result = check_molecule_stability("c1ccccc1")
   print(result["assessment"])       # "Likely stable"
   print(result["stability_score"])  # 100

Checks performed
----------------

1. **Strained small rings** — 3- or 4-membered rings that contain a
   heteroatom (N, O, S …) are flagged.  Aziridines and azetidines are
   significantly more strained than plain cyclopropane / cyclobutane.

2. **Anti-aromatic motifs** — rings whose π-electron count is a multiple
   of 4 (Hückel 4n rule) are thermodynamically destabilised.  Known
   patterns (cyclobutadiene, cyclooctatetraene, pentalene) are matched
   via SMARTS, and π electrons are also counted directly for rings of
   size 4, 8, 12, or 16.

3. **Fused small rings** — two rings of ≤ 4 atoms sharing atoms create
   extreme angle strain (e.g. bicyclo[1.1.0]butane).  These systems can
   be explosive.

4. **Large heterocycles** — rings with ≥ 7 atoms containing heteroatoms
   tend to be conformationally floppy and often unstable.  Very large
   heterocycles (> 10 atoms, ≥ 3 heteroatoms) get an extra penalty.

5. **Carbocations** — positively charged carbon centres are reactive
   intermediates, not isolable species.  The checker distinguishes:

   * *sp2* (``[C+;X3]``) — penalised unless stabilised by an adjacent
     aromatic ring, allylic double bond, or benzylic position.
   * *sp* (``[C+;X2]``) — always penalised heavily.
   * *Primary vs secondary* — primary is worse because fewer alkyl
     groups donate electron density.
   * *Adjacent to EWG* — a carbocation next to F, Cl, Br, I or charged
     N/S/O is the worst case.

6. **Carbenes** — a neutral carbon with two bonds and no hydrogens
   (``[C;X2;H0;+0]``) is extremely reactive.  Extra penalties if the
   carbene sits inside a 3- or 4-membered ring or is adjacent to an
   electron-withdrawing group.

7. **Fused cyclopentane + small hetero ring** — a 5-membered all-carbon
   ring sharing atoms with a 3- or 4-membered hetero ring creates
   significant ring strain.

8. **Physicochemical outliers** — extreme ``logP`` values with
   ``abs(logP) > 10`` or too many rotatable bonds (``> 15``) each
   incur a small penalty.

9. **Aromatic bonus** — aromatic rings stabilise a molecule, so each
   one adds a small bonus (capped at +15 total).

Scoring
-------
After all checks, the score is clamped to 0–100:

* **≥ 80** → ``"Likely stable"``
* **50–79** → ``"Moderately stable"``
* **< 50** → ``"Potentially unstable"``

Entry points
------------

* `check_molecule_stability` — analyse a single SMILES and return a
  0–100 stability score with an issue list.
* `is_valid_smiles` — quick check that a SMILES string parses.

.. automodule:: deepretro.algorithms.stability_checker
   :members:
   :undoc-members:
