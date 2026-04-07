"""r-BRICS balanced molecular splitter for retrosynthesis.

Pipeline
--------
1. **r-BRICS bond discovery** -- ``rBRICS_public.FindrBRICSBonds`` finds every
   bond that the extended BRICS ruleset considers cleavable.  This is a superset
   of standard BRICS (adds fused-ring, long-chain, and heteroatom rules).

2. **Bond classification** -- each bond is tagged as *macrocycle* (in a ring
   with >= ``min_ring_to_cut`` atoms), *small-ring* (protected, never cut), or
   *non-ring*.  The cuttable pool is macrocycle + non-ring bonds.

3. **Combinatorial search with cap** -- we try every way to pick ``n_cuts``
   bonds from the cuttable pool and fragment the molecule.  If the number of
   combinations exceeds ``max_combinations`` (default 50 000), the pool is
   trimmed to the top ``n_cuts + 8`` bonds (macrocycle first, then non-ring
   bonds near heteroatoms) so the search stays fast.

4. **Balance scoring** -- for each candidate split, fragments are scored by
   coefficient of variation (CV = std / mean of heavy-atom counts).  The split
   with the lowest CV wins.  Splits that produce any fragment smaller than
   ``min_size`` are rejected.

5. **Reactive-group replacement** (optional) -- dummy atoms ``[n*]`` on the
   fragments are replaced with chemically meaningful caps (Br, Cl, [H], …)
   looked up from ``reactive_groups.json`` so each fragment is a valid reagent.

Examples
--------
>>> from deepretro.algorithms.splitter import split_molecule
>>> split_molecule("CCCOCc1ccccc1")
['[4*]Cc1ccccc1', '[3*]O[3*]', '[4*]CCC']

>>> from deepretro.algorithms.splitter.rbrics_splitter import replace_dummy_atoms
>>> replace_dummy_atoms(['[4*]CCC', '[3*]O[3*]', '[4*]Cc1ccccc1'])
['BrCCC', 'O', 'BrCc1ccccc1']
"""

import json
import math
import re
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import FragmentOnBonds, GetMolFrags

from .rBRICS_public import FindrBRICSBonds

# Upper bound on combinations before the pool is trimmed.
MAX_COMBINATIONS = 50_000

# How many extra bonds to keep beyond n_cuts when trimming.
POOL_BUFFER = 8

# Isotope offset that signals "use the complement cap" for symmetric cuts.
_COMPLEMENT_OFFSET = 500

_REACTIVE_GROUPS: Optional[dict] = None


def _load_reactive_groups() -> dict:
    """Lazy-load and cache ``reactive_groups.json``."""
    global _REACTIVE_GROUPS
    if _REACTIVE_GROUPS is None:
        p = Path(__file__).with_name("reactive_groups.json")
        with open(p) as fh:
            _REACTIVE_GROUPS = json.load(fh)
    return _REACTIVE_GROUPS


# ------------------------------------------------------------------
# Main API
# ------------------------------------------------------------------

def split_molecule(
    smiles: str,
    n_fragments: Optional[int] = None,
    target_size: int = 15,
    min_size: int = 4,
    min_ring_to_cut: int = 10,
    max_combinations: int = MAX_COMBINATIONS,
) -> List[str]:
    """Split a molecule into balanced fragments for retrosynthesis.

    Parameters
    ----------
    smiles : str
        Input molecule as a SMILES string.
    n_fragments : int, optional
        Target number of fragments.  ``None`` = auto-detect from molecule
        size (roughly ``n_heavy_atoms / target_size``, clamped to 1--7).
    target_size : int, default 15
        Ideal heavy atoms per fragment (used for auto ``n_fragments``).
    min_size : int, default 4
        Reject any split that produces a fragment smaller than this.
    min_ring_to_cut : int, default 10
        Rings smaller than this are protected and never cut.
    max_combinations : int, default 50 000
        If the number of bond subsets exceeds this, the pool is trimmed to
        the top ``n_cuts + 8`` bonds so the search stays fast.  Higher
        values give better results but take longer.

    Returns
    -------
    List[str]
        Fragment SMILES with dummy-atom markers (``[3*]``, ``[4*]``, etc.)
        sorted largest-first by heavy atom count.
        Returns ``[canonical_smiles]`` if the molecule is too small to split.
        Returns ``[]`` if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    canonical = Chem.MolToSmiles(mol)
    n_atoms = mol.GetNumHeavyAtoms()

    if n_fragments is None:
        n_fragments = max(1, min(7, round(n_atoms / target_size)))

    if n_fragments <= 1 or n_atoms <= target_size:
        return [canonical]

    cuttable, _protected = find_cuttable_bonds(mol, min_ring_to_cut)

    if not cuttable:
        return [canonical]

    n_cuts = min(n_fragments - 1, len(cuttable))

    if math.comb(len(cuttable), n_cuts) > max_combinations:
        cuttable = prioritise_pool(cuttable, n_cuts)

    best_cv = float("inf")
    best_frags: Optional[List[str]] = None

    for combo in combinations(range(len(cuttable)), n_cuts):
        bond_indices = [cuttable[i]["bond_idx"] for i in combo]
        dummy_labels = _build_dummy_labels(mol, [cuttable[i] for i in combo])

        try:
            frag_mol = FragmentOnBonds(
                mol, bond_indices, addDummies=True, dummyLabels=dummy_labels
            )
            frag_mols = GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
            sizes = [fm.GetNumHeavyAtoms() for fm in frag_mols]

            if any(s < min_size for s in sizes):
                continue

            cv = calculate_balance_score(sizes)
            if cv < best_cv:
                best_cv = cv
                paired = sorted(
                    ((fm.GetNumHeavyAtoms(), Chem.MolToSmiles(fm)) for fm in frag_mols),
                    reverse=True,
                )
                best_frags = [smi for _, smi in paired]

        except Exception:
            continue

    if best_frags is None:
        return [canonical]

    return best_frags


def replace_dummy_atoms(
    fragments: List[str],
) -> Tuple[List[str], List[dict]]:
    """Replace ``[n*]`` dummy atoms with reactive end-group caps.

    Each dummy carries the BRICS label of the atom it is bonded to (set
    by ``split_molecule`` via ``dummyLabels``).  The cap SMILES is looked
    up from ``reactive_groups.json``.

    For symmetric cuts (e.g. aryl--aryl, 16-16) one dummy is marked
    with an isotope offset (>= 500) so it receives the *complement* cap
    (e.g. boronic acid instead of bromide) to enable Suzuki coupling.

    Parameters
    ----------
    fragments : List[str]
        Fragment SMILES produced by ``split_molecule``.

    Returns
    -------
    capped : List[str]
        Fragment SMILES with dummies replaced by reactive caps.
    info : List[dict]
        Per-dummy replacement log: ``label``, ``cap``, ``complement``,
        ``group``.
    """
    lookup = _load_reactive_groups()
    replacements = lookup["replacements"]
    complements = lookup.get("complements", {})

    capped: List[str] = []
    info: List[dict] = []

    dummy_re = re.compile(r"\[(\d+)\*\]")

    for frag_smi in fragments:
        new_smi = frag_smi

        for m in dummy_re.finditer(frag_smi):
            isotope = int(m.group(1))
            use_complement = isotope >= _COMPLEMENT_OFFSET
            label = str(isotope - _COMPLEMENT_OFFSET) if use_complement else str(isotope)

            entry = replacements.get(label)
            if entry is None or not entry["cap"]:
                continue

            default_cap = entry["cap"]

            if use_complement:
                comp = complements.get(default_cap)
                cap = comp["cap"] if comp else default_cap
            else:
                cap = default_cap

            new_smi = new_smi.replace(m.group(0), cap, 1)
            info.append({
                "label": label,
                "cap": cap,
                "complement": use_complement,
                "group": entry["group"],
            })

        mol = Chem.MolFromSmiles(new_smi)
        if mol is not None:
            capped.append(Chem.MolToSmiles(mol))
        else:
            capped.append(frag_smi)

    return capped, info


# ------------------------------------------------------------------
# Bond discovery & helpers
# ------------------------------------------------------------------

def find_cuttable_bonds(
    mol: Chem.Mol,
    min_ring_to_cut: int = 10,
) -> Tuple[List[Dict], List[Dict]]:
    """Find r-BRICS bonds and classify by ring membership.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    min_ring_to_cut : int, default 10
        Minimum ring size to allow cutting.

    Returns
    -------
    cuttable : List[Dict]
        Macrocycle + non-ring bonds.  Each dict has keys
        ``bond_idx``, ``atoms``, ``labels``, ``type``, ``symbols``.
    protected : List[Dict]
        Bonds in small rings that will not be cut.
    """
    raw_bonds = list(FindrBRICSBonds(mol))
    atom_rings = list(mol.GetRingInfo().AtomRings())

    cuttable: List[Dict] = []
    protected: List[Dict] = []
    seen: set = set()

    for (a1, a2), (t1, t2) in raw_bonds:
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is None:
            continue
        bid = bond.GetIdx()
        if bid in seen:
            continue
        seen.add(bid)

        info: Dict = {
            "bond_idx": bid,
            "atoms": (a1, a2),
            "labels": (t1, t2),
            "type": f"L{t1}-L{t2}",
            "symbols": (
                f"{mol.GetAtomWithIdx(a1).GetSymbol()}({a1})-"
                f"{mol.GetAtomWithIdx(a2).GetSymbol()}({a2})"
            ),
        }

        if not bond.IsInRing():
            cuttable.append(info)
            continue

        in_small = False
        in_macro = False
        for ring_atoms in atom_rings:
            if a1 in ring_atoms and a2 in ring_atoms:
                if len(ring_atoms) < min_ring_to_cut:
                    in_small = True
                else:
                    in_macro = True

        if in_small:
            protected.append(info)
        elif in_macro:
            cuttable.insert(0, info)
        else:
            protected.append(info)

    return cuttable, protected


def calculate_balance_score(sizes: List[int]) -> float:
    """Coefficient of variation of fragment sizes.  Lower = more balanced.

    Parameters
    ----------
    sizes : List[int]
        Heavy-atom counts per fragment.

    Returns
    -------
    float
        ``std / mean``.  Returns ``inf`` for < 2 fragments or zero mean.

    Examples
    --------
    >>> calculate_balance_score([10, 10, 10])
    0.0
    >>> round(calculate_balance_score([10, 20]), 4)
    0.3333
    """
    if len(sizes) < 2:
        return float("inf")
    mean = sum(sizes) / len(sizes)
    if mean == 0:
        return float("inf")
    variance = sum((s - mean) ** 2 for s in sizes) / len(sizes)
    return math.sqrt(variance) / mean


def prioritise_pool(cuttable: List[Dict], n_cuts: int) -> List[Dict]:
    """Trim the bond pool when combinations would exceed the safety cap.

    Macrocycle bonds (already at the front of *cuttable* from
    ``find_cuttable_bonds``) are kept first.  Among non-ring bonds,
    those adjacent to O or N are prioritised because they tend to
    correspond to chemically meaningful disconnections (ethers, amines,
    amides, esters).  The pool is then truncated to ``n_cuts + 8`` so
    the combinatorial search stays under ~tens-of-thousands of subsets.
    """
    cuttable = sorted(
        cuttable,
        key=lambda b: 0 if "O" in b["symbols"] or "N" in b["symbols"] else 1,
    )
    return cuttable[: n_cuts + POOL_BUFFER]


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _build_dummy_labels(
    mol: Chem.Mol,
    bonds: List[Dict],
) -> List[Tuple[int, int]]:
    """Build ``dummyLabels`` for ``FragmentOnBonds``.

    Each pair ``(label_for_begin_atom, label_for_end_atom)`` uses the
    BRICS label of the atom that the dummy is *bonded to* (same
    convention as ``BreakrBRICSBonds``).  For symmetric cuts (both
    labels equal), the end-atom side is offset by 500 so
    ``replace_dummy_atoms`` can apply the complement cap.
    """
    labels: List[Tuple[int, int]] = []

    for b in bonds:
        a1, a2 = b["atoms"]
        t1, t2 = b["labels"]
        bond = mol.GetBondBetweenAtoms(a1, a2)
        begin = bond.GetBeginAtomIdx()

        is_symmetric = t1 == t2

        if begin == a1:
            lab_begin = int(t1)
            lab_end = int(t2) + (_COMPLEMENT_OFFSET if is_symmetric else 0)
        else:
            lab_begin = int(t2) + (_COMPLEMENT_OFFSET if is_symmetric else 0)
            lab_end = int(t1)

        labels.append((lab_begin, lab_end))

    return labels
