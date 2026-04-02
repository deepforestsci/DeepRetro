"""Heuristic stability checker for molecules in retrosynthetic pathways.

When an LLM proposes reactant molecules, some of them may be chemically
unstable — strained small rings, anti-aromatic systems, reactive
intermediates like carbocations or carbenes, etc.  This module catches
those problems by inspecting the molecular graph with RDKit descriptors
and SMARTS pattern matching.  No trained model is required.

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
"""

from __future__ import annotations

from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import (
    CalcNumAliphaticCarbocycles,
    CalcNumAliphaticHeterocycles,
    CalcNumAliphaticRings,
    CalcNumAromaticCarbocycles,
    CalcNumAromaticHeterocycles,
    CalcNumAromaticRings,
    CalcNumBridgeheadAtoms,
)

# Helpers


def is_valid_smiles(smiles: str) -> bool:
    """Check whether a SMILES string can be parsed by RDKit.

    Parameters
    ----------
    smiles : str
        SMILES string to validate.

    Returns
    -------
    bool
        ``True`` if RDKit can build a molecule from *smiles*.

    Examples
    --------
    >>> is_valid_smiles("CCO")
    True
    >>> is_valid_smiles("not_a_molecule")
    False
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return False
    return mol is not None


# Anti-aromatic / reactive-intermediate SMARTS (compiled once)

_ANTI_AROMATIC_PATTERNS: list[tuple[Chem.Mol | None, str]] = [
    (Chem.MolFromSmarts("c1ccc1"), "cyclobutadiene-like (anti-aromatic)"),
    (
        Chem.MolFromSmarts("C1=CC=CC=CC=C1"),
        "cyclooctatetraene-like (potential anti-aromatic)",
    ),
    (Chem.MolFromSmarts("c1cc2cccc2c1"), "pentalene-like (potential anti-aromatic)"),
]

_SP2_CARBOCATION = Chem.MolFromSmarts("[C+;X3]")
_SP_CARBOCATION = Chem.MolFromSmarts("[C+;X2]")
_UNSTABILISED_CARBOCATION = Chem.MolFromSmarts("[C+][F,Cl,Br,I,N+,S+,O+]")
_PRIMARY_CARBOCATION = Chem.MolFromSmarts("[CH2+][#6]")
_SECONDARY_CARBOCATION = Chem.MolFromSmarts("[CH+]([#6])[#6]")
_ALLYLIC_CARBOCATION = Chem.MolFromSmarts("[C+]-[C]=[C]")
_BENZYLIC_CARBOCATION = Chem.MolFromSmarts("[C+]-c")

_SINGLET_CARBENE = Chem.MolFromSmarts("[C;X2;H0;+0]")
_UNSTABLE_CARBENE = Chem.MolFromSmarts("[C;X2;H0;+0][F,Cl,Br,I,N+,S+,O+]")
_RING3_CARBENE = Chem.MolFromSmarts("[C;X2;H0;+0]1[C,N,O][C,N,O]1")
_RING4_CARBENE = Chem.MolFromSmarts("[C;X2;H0;+0]1[C,N,O][C,N,O][C,N,O]1")

_CYCLOPENTANE = Chem.MolFromSmarts("C1CCCC1")

# Core


def check_molecule_stability(smiles: str) -> dict[str, Any]:
    """Assess the stability of a molecule from its SMILES string.

    Parses the molecule with RDKit and runs nine heuristic checks
    (see the module docstring for a full description of each one):
    strained small rings, anti-aromatic motifs, fused small rings,
    large heterocycles, carbocations, carbenes, fused cyclopentane
    systems, physicochemical outliers, and an aromatic-ring bonus.
    Each problem subtracts from a base score of 100; the final score
    is clamped to 0–100.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule to assess.

    Returns
    -------
    dict[str, Any]
        Keys returned:

        * ``valid_structure`` (bool) — whether RDKit could parse the
          SMILES at all.
        * ``stability_score`` (int, 0–100) — overall stability rating.
        * ``issues`` (list[str]) — plain-English descriptions of every
          problem found (e.g. ``"Three-membered heterocycle
          (potentially unstable)"``).
        * ``metrics`` (dict) — molecular weight, logP, H-bond donors /
          acceptors, rotatable bonds.
        * ``ring_data`` (dict) — ring counts broken down by type
          (aliphatic / aromatic, carbocycle / heterocycle, bridgehead
          atoms, etc.).
        * ``atom_data`` (dict) — total atoms, bonds, heavy atoms,
          aromatic vs aliphatic counts.
        * ``assessment`` (str) — one of ``"Likely stable"``,
          ``"Moderately stable"``, or ``"Potentially unstable"``.

    Examples
    --------
    >>> res = check_molecule_stability("c1ccccc1")  # benzene
    >>> res["assessment"]
    'Likely stable'
    >>> res = check_molecule_stability("[CH2+]C")    # ethyl cation
    >>> res["assessment"]
    'Potentially unstable'
    """
    results: dict[str, Any] = {
        "valid_structure": False,
        "stability_score": 0,
        "issues": [],
        "metrics": {},
        "ring_data": {},
        "atom_data": {},
        "assessment": "",
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        results["issues"].append("Invalid SMILES string or cannot be parsed")
        return results

    results["valid_structure"] = True

    # basic physicochemical properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    results["metrics"] = {
        "molecular_weight": mw,
        "logP": logp,
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": rotatable_bonds,
    }

    # ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()

    num_aromatic_atoms = len(mol.GetAromaticAtoms())
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    num_aliphatic_carbocycles = CalcNumAliphaticCarbocycles(mol)
    num_aliphatic_heterocycles = CalcNumAliphaticHeterocycles(mol)
    num_aliphatic_rings = CalcNumAliphaticRings(mol)
    num_aromatic_carbocycles = CalcNumAromaticCarbocycles(mol)
    num_aromatic_heterocycles = CalcNumAromaticHeterocycles(mol)
    num_aromatic_rings = CalcNumAromaticRings(mol)
    num_bridgehead_atoms = CalcNumBridgeheadAtoms(mol)

    results["ring_data"] = {
        "num_rings": len(atom_rings),
        "atom_rings": [list(r) for r in atom_rings],
        "bond_rings": [list(r) for r in bond_rings],
        "num_aliphatic_carbocycles": num_aliphatic_carbocycles,
        "num_aliphatic_heterocycles": num_aliphatic_heterocycles,
        "num_aliphatic_rings": num_aliphatic_rings,
        "num_aromatic_carbocycles": num_aromatic_carbocycles,
        "num_aromatic_heterocycles": num_aromatic_heterocycles,
        "num_aromatic_rings": num_aromatic_rings,
        "num_bridgehead_atoms": num_bridgehead_atoms,
    }

    results["atom_data"] = {
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "num_heavy_atoms": num_heavy_atoms,
        "num_heavy_bonds": mol.GetNumBonds(),
        "num_aromatic_atoms": num_aromatic_atoms,
        "num_aliphatic_atoms": num_heavy_atoms - num_aromatic_atoms,
    }

    # strained small rings
    for ring in atom_rings:
        if len(ring) < 3:
            results["issues"].append(f"Highly strained ring of size {len(ring)}")
        elif len(ring) == 3 and any(
            mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in ring
        ):
            results["issues"].append(
                "Three-membered heterocycle (potentially unstable)"
            )
        elif len(ring) == 4 and any(
            mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in ring
        ):
            results["issues"].append("Four-membered heterocycle (potentially unstable)")

    # anti-aromatic motifs
    for patt, name in _ANTI_AROMATIC_PATTERNS:
        if patt and mol.HasSubstructMatch(patt):
            results["issues"].append(f"Contains {name} motif")

    for ring in atom_rings:
        if len(ring) in (4, 8, 12, 16):
            double_bond_count = 0.0
            ring_atoms_obj = [mol.GetAtomWithIdx(i) for i in ring]
            for atom in ring_atoms_obj:
                if atom.GetIsAromatic():
                    double_bond_count += 0.5
                for bond in atom.GetBonds():
                    other = bond.GetOtherAtomIdx(atom.GetIdx())
                    if other in ring and bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bond_count += 0.5
            pi_electrons = int(double_bond_count * 2)
            if pi_electrons > 0 and pi_electrons % 4 == 0:
                results["issues"].append(
                    f"{len(ring)}-membered ring with {pi_electrons} "
                    f"\u03c0 electrons (potential anti-aromatic)"
                )

    # fused small rings
    small_rings = [set(r) for r in atom_rings if len(r) <= 4]
    fused_small_rings_detected = False
    for i in range(len(small_rings)):
        for j in range(i + 1, len(small_rings)):
            if small_rings[i] & small_rings[j]:
                fused_small_rings_detected = True
                results["issues"].append(
                    f"Fused {len(small_rings[i])} and "
                    f"{len(small_rings[j])}-membered rings "
                    f"(highly strained system)"
                )

    if fused_small_rings_detected:
        results["issues"].append(
            "WARNING: Fused small rings create highly strained "
            "and potentially explosive compounds"
        )

    # large heterocycles
    for ring in atom_rings:
        if len(ring) >= 7:
            ring_atoms_obj = [mol.GetAtomWithIdx(i) for i in ring]
            heteroatoms = [a for a in ring_atoms_obj if a.GetSymbol() not in ("C", "H")]
            if heteroatoms:
                symbols = {a.GetSymbol() for a in heteroatoms}
                results["issues"].append(
                    f"Large ({len(ring)}-membered) heterocycle with "
                    f"{', '.join(symbols)} (potentially unstable)"
                )
                if len(ring) > 10 and len(heteroatoms) >= 3:
                    results["issues"].append(
                        "Very large heterocycle with multiple heteroatoms "
                        "(significant stability concern)"
                    )

    # bridgehead strain
    if num_bridgehead_atoms > 0:
        if any(len(r) <= 4 for r in atom_rings) and num_bridgehead_atoms >= 2:
            results["issues"].append("Strained polycyclic system with bridgehead atoms")

    # scoring (base 100, subtract penalties)
    score = 100
    score -= len(results["issues"]) * 15

    # carbocations
    score = check_carbocations(mol, results, score)

    # carbenes
    score = check_carbenes(mol, results, score)

    # fused cyclopentane + small hetero ring
    score = check_fused_cyclopentane(mol, atom_rings, results, score)

    # physicochemical penalties
    if abs(logp) > 10:
        score -= 10
    if rotatable_bonds > 15:
        score -= 10

    # aromatic bonus
    if num_aromatic_rings > 0:
        score += min(num_aromatic_rings * 5, 15)

    # aliphatic heterocycle penalty
    if num_aliphatic_heterocycles > 0 and num_aliphatic_rings <= 3:
        score -= 5 * num_aliphatic_heterocycles

    # bridgehead penalty
    if num_bridgehead_atoms > 0 and any(len(r) <= 5 for r in atom_rings):
        score -= num_bridgehead_atoms * 5

    if fused_small_rings_detected:
        score -= 30

    # anti-aromatic penalty
    for issue in results["issues"]:
        if "anti-aromatic" in issue:
            score -= 25
        elif "\u03c0 electrons" in issue:
            score -= 20

    # large heterocycle penalty
    large_hetero = sum(
        1
        for iss in results["issues"]
        if "heterocycle" in iss and "large" in iss.lower()
    )
    if large_hetero > 0:
        score -= large_hetero * 10

    # clamp and label
    score = max(0, min(100, score))
    results["stability_score"] = score

    if score >= 80:
        results["assessment"] = "Likely stable"
    elif score >= 50:
        results["assessment"] = "Moderately stable"
    else:
        results["assessment"] = "Potentially unstable"

    return results


def check_carbocations(
    mol: Chem.Mol,
    results: dict[str, Any],
    score: int,
) -> int:
    """Detect carbocation intermediates and apply score penalties.

    Carbocations are positively charged carbon centres i.e. reactive
    intermediates that cannot be bottled.  The function uses SMARTS
    pattern matching to find them and then checks whether the charge
    is stabilised by resonance (allylic or benzylic position) or by
    neighbouring aromatic atoms.  Unstabilised and primary carbocations
    get the heaviest penalties; stabilised ones get a small bonus.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object already parsed from SMILES.
    results : dict[str, Any]
        Accumulator — detected issues are appended to ``results["issues"]``.
    score : int
        Running stability score to adjust.

    Returns
    -------
    int
        Updated stability score after carbocation penalties / bonuses.

    Examples
    --------
    >>> from rdkit import Chem
    >>> from deepretro.algorithms import check_carbocations
    >>> mol = Chem.MolFromSmiles("[CH2+]C")  # ethyl cation
    >>> results = {"issues": []}
    >>> new_score = check_carbocations(mol, results, 100)
    >>> "Contains primary carbocation (highly unstable)" in results["issues"]
    True
    >>> new_score < 100
    True
    """
    if _SP2_CARBOCATION and mol.HasSubstructMatch(_SP2_CARBOCATION):
        for match in mol.GetSubstructMatches(_SP2_CARBOCATION):
            carbon = mol.GetAtomWithIdx(match[0])
            neighbours = [mol.GetAtomWithIdx(n.GetIdx()) for n in carbon.GetNeighbors()]
            is_allylic = _ALLYLIC_CARBOCATION and mol.HasSubstructMatch(
                _ALLYLIC_CARBOCATION
            )
            is_benzylic = _BENZYLIC_CARBOCATION and mol.HasSubstructMatch(
                _BENZYLIC_CARBOCATION
            )
            if any(a.GetIsAromatic() for a in neighbours) or is_allylic or is_benzylic:
                score += 10
            else:
                results["issues"].append(
                    "Contains non-stabilized sp2 carbocation "
                    "(highly unstable intermediate)"
                )
                score -= 25

    if _SP_CARBOCATION and mol.HasSubstructMatch(_SP_CARBOCATION):
        results["issues"].append(
            "Contains sp carbocation (highly unstable intermediate)"
        )
        score -= 30

    if _UNSTABILISED_CARBOCATION and mol.HasSubstructMatch(_UNSTABILISED_CARBOCATION):
        results["issues"].append(
            "Contains carbocation adjacent to electron-withdrawing "
            "group (extremely unstable)"
        )
        score -= 35

    if _PRIMARY_CARBOCATION and mol.HasSubstructMatch(_PRIMARY_CARBOCATION):
        results["issues"].append("Contains primary carbocation (highly unstable)")
        score -= 30
    elif _SECONDARY_CARBOCATION and mol.HasSubstructMatch(_SECONDARY_CARBOCATION):
        results["issues"].append(
            "Contains secondary carbocation (unstable intermediate)"
        )
        score -= 25

    return score


def check_carbenes(
    mol: Chem.Mol,
    results: dict[str, Any],
    score: int,
) -> int:
    """Detect carbene intermediates and apply score penalties.

    A carbene is a neutral carbon with only two bonds and no hydrogens
    i.e. an extremely reactive species that usually exists only as a
    fleeting intermediate.  Additional penalties stack if the carbene
    is inside a strained 3- or 4-membered ring, or sits next to an
    electron-withdrawing group (halogens, charged heteroatoms).

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object already parsed from SMILES.
    results : dict[str, Any]
        Accumulator — detected issues are appended to ``results["issues"]``.
    score : int
        Running stability score to adjust.

    Returns
    -------
    int
        Updated stability score after carbene penalties.

    Examples
    --------
    >>> from rdkit import Chem
    >>> from deepretro.algorithms import check_carbenes
    >>> mol = Chem.MolFromSmiles("[C]1CC1")  # carbene in 3-membered ring
    >>> results = {"issues": []}
    >>> new_score = check_carbenes(mol, results, 100)
    >>> any("carbene" in i for i in results["issues"])
    True
    >>> new_score < 100
    True
    """
    if _SINGLET_CARBENE and mol.HasSubstructMatch(_SINGLET_CARBENE):
        results["issues"].append("Contains carbene (highly reactive intermediate)")
        score -= 35

        if _UNSTABLE_CARBENE and mol.HasSubstructMatch(_UNSTABLE_CARBENE):
            results["issues"].append(
                "Contains carbene adjacent to electron-withdrawing "
                "group (extremely unstable)"
            )
            score -= 10

        if _RING3_CARBENE and mol.HasSubstructMatch(_RING3_CARBENE):
            results["issues"].append(
                "Contains carbene in 3-membered ring (extremely unstable)"
            )
            score -= 15

        if _RING4_CARBENE and mol.HasSubstructMatch(_RING4_CARBENE):
            results["issues"].append(
                "Contains carbene in 4-membered ring (highly unstable)"
            )
            score -= 10

    return score


def check_fused_cyclopentane(
    mol: Chem.Mol,
    atom_rings: tuple,
    results: dict[str, Any],
    score: int,
) -> int:
    """Detect 5-membered carbon rings fused with small hetero rings.

    A cyclopentane ring sharing atoms with a 3- or 4-membered ring
    that contains a heteroatom (N, O, S …) creates significant angle
    strain, for example 1,2-epoxycyclopentane.  Each such fusion
    incurs a heavy penalty (-40).

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object already parsed from SMILES.
    atom_rings : tuple
        Ring atom-index tuples from ``mol.GetRingInfo().AtomRings()``.
    results : dict[str, Any]
        Accumulator — detected issues are appended to ``results["issues"]``.
    score : int
        Running stability score to adjust.

    Returns
    -------
    int
        Updated stability score after fused-ring penalties.

    Examples
    --------
    >>> from rdkit import Chem
    >>> from deepretro.algorithms import check_fused_cyclopentane
    >>> mol = Chem.MolFromSmiles("C1CC2OCC12")  # epoxycyclopentane
    >>> rings = mol.GetRingInfo().AtomRings()
    >>> results = {"issues": []}
    >>> new_score = check_fused_cyclopentane(mol, rings, results, 100)
    >>> any("strained system" in i for i in results["issues"])
    True
    >>> new_score < 100
    True
    """
    if not (_CYCLOPENTANE and mol.HasSubstructMatch(_CYCLOPENTANE)):
        return score

    five_mem = [
        r
        for r in atom_rings
        if len(r) == 5 and all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)
    ]
    small_hetero = [
        r
        for r in atom_rings
        if len(r) in (3, 4) and any(mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in r)
    ]

    for fr in five_mem:
        fr_set = set(fr)
        for sr in small_hetero:
            if fr_set & set(sr):
                hetero = {
                    mol.GetAtomWithIdx(i).GetSymbol()
                    for i in sr
                    if mol.GetAtomWithIdx(i).GetSymbol() != "C"
                }
                results["issues"].append(
                    f"5-membered carbon ring fused with "
                    f"{len(sr)}-membered ring containing "
                    f"{', '.join(hetero)} (strained system)"
                )
                score -= 40

    return score
