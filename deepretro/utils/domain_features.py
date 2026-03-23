"""Domain feature extraction utilities for reaction-step featurization.

Extended feature set (v2) — adds ring-topology, functional-group,
similarity, and structural-change features motivated by chemist
feedback on common LLM hallucination modes:

* H1  Ring size changes  (6→5, 6→7, 7→4)
* H2  Atom-count drift / chain-length changes
* H3  Functional-group swaps  (CH₃↔OH, CH₃↔OCH₃)
* H4  Reduction loops  (molecule barely changes between steps)
* H8  Unstable intermediates  (3/4-membered rings)
* H9  Sudden size collapse  (molecule shrinks drastically)
"""

import numpy as np
from collections import Counter
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors

# ---------------------------------------------------------------------------
# Total count of hand-crafted features — must stay in sync with the
# np.array([...]) at the end of extract_domain_features_single().
# Original 15  +  12 new  =  27
# ---------------------------------------------------------------------------
NUM_DOMAIN_FEATURES = 27

# ---------------------------------------------------------------------------
# SMARTS patterns for functional-group counting  (H3 detection)
# ---------------------------------------------------------------------------
_FG_SMARTS: dict[str, str] = {
    "hydroxyl":  "[OX2H]",           # -OH
    "carbonyl":  "[CX3]=[OX1]",      # C=O  (ketone / aldehyde / acid / ester)
    "methyl":    "[CH3;!$([CH3]~[!#6])]",  # -CH3 bonded only to C
    "amine":     "[NX3;H2,H1;!$(NC=O)]",   # primary / secondary amine
}
_FG_PATTERNS: dict[str, Chem.Mol] = {}


def _get_fg_pattern(name: str) -> Chem.Mol | None:
    """Lazy-compile SMARTS only once."""
    if name not in _FG_PATTERNS:
        _FG_PATTERNS[name] = Chem.MolFromSmarts(_FG_SMARTS[name])
    return _FG_PATTERNS[name]


def _count_fg(mol: Chem.Mol, name: str) -> int:
    """Count non-overlapping matches of a named functional-group SMARTS."""
    pat = _get_fg_pattern(name)
    if pat is None:
        return 0
    return len(mol.GetSubstructMatches(pat))


# ---------------------------------------------------------------------------
# Ring helpers  (H1 / H8 detection)
# ---------------------------------------------------------------------------

def _ring_sizes(mol: Chem.Mol) -> list[int]:
    """Return sorted list of ring sizes from SSSR."""
    return sorted(len(r) for r in Chem.GetSSSR(mol))


def _ring_size_distance(sizes_a: list[int], sizes_b: list[int]) -> float:
    """Histogram-style distance between two ring-size distributions.

    Builds a count vector for ring sizes 3..12 for each side and
    returns the L1 (Manhattan) distance.  This captures not just the
    *number* of rings changing, but *which* sizes changed.
    """
    def _hist(sizes):
        h = np.zeros(10, dtype=float)  # indices 0..9 → ring sizes 3..12
        for s in sizes:
            idx = min(max(s - 3, 0), 9)
            h[idx] += 1
        return h
    return float(np.sum(np.abs(_hist(sizes_a) - _hist(sizes_b))))


# ---------------------------------------------------------------------------
# Tanimoto similarity  (H4 / H9 detection)
# ---------------------------------------------------------------------------

def _tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol,
              radius: int = 2, n_bits: int = 2048) -> float:
    """Morgan-fingerprint Tanimoto between two molecules."""
    fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius, nBits=n_bits)
    fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius, nBits=n_bits)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


# ---------------------------------------------------------------------------
# Main feature extractor
# ---------------------------------------------------------------------------

def extract_domain_features_single(
    product_smiles: str,
    reactants_smiles: str,
) -> np.ndarray:
    """Extract hand-crafted domain features for one product → reactant pair.

    **Original 15 features** (indices 0-14):
        Atom-count deltas (C, N, O, Cl, Br), bond / ring / aromaticity
        deltas, molecular-weight deltas, num-reactants, absolute counts.

    **New 12 features** (indices 15-26), targeting chemist-identified
    hallucination modes:

    15. ``max_ring_delta``      – max ring size in reactant minus product
    16. ``min_ring_delta``      – min ring size in reactant minus product
    17. ``small_ring_count_r``  – number of 3-or-4-membered rings in reactant
    18. ``ring_hist_distance``  – L1 distance of ring-size histograms
    19. ``tanimoto_sim``        – Morgan-FP Tanimoto (product vs combined reactant)
    20. ``heavy_atom_ratio``    – reactant heavy atoms / product heavy atoms
    21. ``hydroxyl_delta``      – reactant -OH count minus product -OH count
    22. ``carbonyl_delta``      – reactant C=O count minus product C=O count
    23. ``methyl_delta``        – reactant -CH₃ count minus product -CH₃ count
    24. ``amine_delta``         – reactant amine count minus product amine count
    25. ``stereocenter_delta``  – reactant chiral centres minus product
    26. ``unsaturation_delta``  – reactant degree-of-unsat. minus product

    Parameters
    ----------
    product_smiles : str
        SMILES of the target product.
    reactants_smiles : str
        SMILES of the proposed reactants (dot-separated when multiple).

    Returns
    -------
    features : np.ndarray, shape (NUM_DOMAIN_FEATURES,)
        1-D feature vector.  Returns a NaN vector on any parsing failure.

    Examples
    --------
    >>> from deepretro.utils import extract_domain_features_single
    >>> feats = extract_domain_features_single("CCO", "CC.O")
    >>> feats.shape
    (27,)
    """
    try:
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol is None:
            return np.full(NUM_DOMAIN_FEATURES, np.nan)

        # Reactants may be multiple molecules joined by '.'
        reactant_mols = [Chem.MolFromSmiles(r)
                         for r in reactants_smiles.split(".")]

        # ── Product-side descriptors ──────────────────────────────────
        p_atoms = Counter(a.GetSymbol() for a in product_mol.GetAtoms())
        p_bonds = product_mol.GetNumBonds()
        p_rings = len(Chem.GetSSSR(product_mol))
        p_arom  = sum(1 for a in product_mol.GetAtoms() if a.GetIsAromatic())
        p_mw    = Descriptors.MolWt(product_mol)

        p_ring_sizes   = _ring_sizes(product_mol)
        p_heavy        = product_mol.GetNumHeavyAtoms()
        p_stereo       = len(Chem.FindMolChiralCenters(product_mol,
                                                        includeUnassigned=True))
        p_unsat        = rdMolDescriptors.CalcNumHBA(product_mol)  # proxy via HBA
        # A better proxy: classical degree of unsaturation = (2C+2+N-H) / 2
        # but H count is implicit in RDKit; use ring+double-bond count instead
        p_dbu = p_rings + sum(
            1 for b in product_mol.GetBonds()
            if b.GetBondTypeAsDouble() >= 2.0
        )

        # ── Reactant-side descriptors (aggregated) ────────────────────
        r_atoms: Counter   = Counter()
        r_bonds            = 0
        r_rings            = 0
        r_arom             = 0
        r_mw               = 0.0
        n_valid            = 0
        r_ring_sizes_all: list[int] = []
        r_heavy            = 0
        r_stereo           = 0
        r_dbu              = 0

        for mol in reactant_mols:
            if mol is None:
                continue
            r_atoms += Counter(a.GetSymbol() for a in mol.GetAtoms())
            r_bonds += mol.GetNumBonds()
            mol_rings = Chem.GetSSSR(mol)
            r_rings += len(mol_rings)
            r_arom  += sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
            r_mw    += Descriptors.MolWt(mol)
            n_valid += 1
            r_ring_sizes_all.extend(sorted(len(r) for r in mol_rings))
            r_heavy  += mol.GetNumHeavyAtoms()
            r_stereo += len(Chem.FindMolChiralCenters(mol,
                                                       includeUnassigned=True))
            r_dbu += len(mol_rings) + sum(
                1 for b in mol.GetBonds()
                if b.GetBondTypeAsDouble() >= 2.0
            )

        # Combined reactant mol for Tanimoto (parse the whole string as one)
        combined_reactant_mol = Chem.MolFromSmiles(reactants_smiles)

        # ── NEW: ring topology features (H1, H8) ─────────────────────
        r_ring_sizes_all_sorted = sorted(r_ring_sizes_all)
        max_ring_delta = (
            (max(r_ring_sizes_all_sorted) if r_ring_sizes_all_sorted else 0)
            - (max(p_ring_sizes) if p_ring_sizes else 0)
        )
        min_ring_delta = (
            (min(r_ring_sizes_all_sorted) if r_ring_sizes_all_sorted else 0)
            - (min(p_ring_sizes) if p_ring_sizes else 0)
        )
        small_ring_count_r = sum(1 for s in r_ring_sizes_all_sorted if s <= 4)
        ring_hist_dist = _ring_size_distance(r_ring_sizes_all_sorted,
                                              p_ring_sizes)

        # ── NEW: similarity / size features (H4, H9) ─────────────────
        if combined_reactant_mol is not None:
            tanimoto_sim = _tanimoto(product_mol, combined_reactant_mol)
        else:
            tanimoto_sim = 0.0

        heavy_atom_ratio = r_heavy / max(p_heavy, 1)

        # ── NEW: functional-group deltas (H3) ────────────────────────
        # Aggregate FG counts across all valid reactant molecules
        r_fg: dict[str, int] = {fg: 0 for fg in _FG_SMARTS}
        for mol in reactant_mols:
            if mol is not None:
                for fg in _FG_SMARTS:
                    r_fg[fg] += _count_fg(mol, fg)

        p_fg: dict[str, int] = {fg: _count_fg(product_mol, fg)
                                 for fg in _FG_SMARTS}

        hydroxyl_delta = r_fg["hydroxyl"] - p_fg["hydroxyl"]
        carbonyl_delta = r_fg["carbonyl"] - p_fg["carbonyl"]
        methyl_delta   = r_fg["methyl"]   - p_fg["methyl"]
        amine_delta    = r_fg["amine"]    - p_fg["amine"]

        # ── NEW: stereocenter & unsaturation deltas (H3, H4) ─────────
        stereocenter_delta = r_stereo - p_stereo
        unsaturation_delta = r_dbu - p_dbu

        # ── Assemble full feature vector ──────────────────────────────
        return np.array([
            # --- original 15 (indices 0-14) ---
            float(r_atoms.get('C', 0) - p_atoms.get('C', 0)),    #  0: C delta
            float(r_atoms.get('N', 0) - p_atoms.get('N', 0)),    #  1: N delta
            float(r_atoms.get('O', 0) - p_atoms.get('O', 0)),    #  2: O delta
            float(r_atoms.get('Cl', 0) - p_atoms.get('Cl', 0)),  #  3: Cl delta
            float(r_atoms.get('Br', 0) - p_atoms.get('Br', 0)),  #  4: Br delta
            float(r_bonds - p_bonds),                              #  5: bond count delta
            float(r_rings - p_rings),                              #  6: ring count delta
            float(r_arom  - p_arom),                               #  7: aromatic atom delta
            float(r_mw    - p_mw),                                 #  8: MW delta
            float(n_valid),                                        #  9: num valid reactants
            float(sum(r_atoms.values()) - sum(p_atoms.values())),  # 10: total heavy-atom delta
            float(p_mw),                                           # 11: product MW
            float(r_mw),                                           # 12: reactant MW (sum)
            float(p_rings),                                        # 13: product ring count
            float(r_rings),                                        # 14: reactant ring count
            # --- NEW 12 features (indices 15-26) ---
            float(max_ring_delta),       # 15: max ring size delta  (H1)
            float(min_ring_delta),       # 16: min ring size delta  (H1)
            float(small_ring_count_r),   # 17: #(3/4-membered rings) in reactant (H8)
            float(ring_hist_dist),       # 18: ring-size histogram L1 distance (H1)
            float(tanimoto_sim),         # 19: Morgan Tanimoto similarity (H4/H9)
            float(heavy_atom_ratio),     # 20: reactant/product heavy-atom ratio (H9)
            float(hydroxyl_delta),       # 21: -OH count delta (H3)
            float(carbonyl_delta),       # 22: C=O count delta (H3)
            float(methyl_delta),         # 23: -CH3 count delta (H3)
            float(amine_delta),          # 24: amine count delta (H3)
            float(stereocenter_delta),   # 25: chiral-centre delta (H3)
            float(unsaturation_delta),   # 26: degree-of-unsaturation delta (H4)
        ])
    except Exception:
        return np.full(NUM_DOMAIN_FEATURES, np.nan)
