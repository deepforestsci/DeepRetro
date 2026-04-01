"""Domain feature extraction utilities for reaction-step featurization."""

import numpy as np
from collections import Counter
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors

NUM_DOMAIN_FEATURES = 27

_FG_SMARTS: dict[str, str] = {
    "hydroxyl":  "[OX2H]",
    "carbonyl":  "[CX3]=[OX1]",
    "methyl":    "[CH3;!$([CH3]~[!#6])]",
    "amine":     "[NX3;H2,H1;!$(NC=O)]",
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


def _ring_sizes(mol: Chem.Mol) -> list[int]:
    """Return sorted list of ring sizes from SSSR."""
    return sorted(len(r) for r in Chem.GetSSSR(mol))


def _ring_size_distance(sizes_a: list[int], sizes_b: list[int]) -> float:
    """Histogram-style L1 distance between two ring-size distributions."""
    def _hist(sizes):
        h = np.zeros(10, dtype=float)
        for s in sizes:
            idx = min(max(s - 3, 0), 9)
            h[idx] += 1
        return h
    return float(np.sum(np.abs(_hist(sizes_a) - _hist(sizes_b))))


def _tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol,
              radius: int = 2, n_bits: int = 2048) -> float:
    """Morgan-fingerprint Tanimoto between two molecules."""
    fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius, nBits=n_bits)
    fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius, nBits=n_bits)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


def extract_domain_features_single(
    product_smiles: str,
    reactants_smiles: str,
) -> np.ndarray:
    """Extract domain features for one product → reactant pair.

    Computes a 27-dimensional vector covering atom-count deltas
    (C, N, O, Cl, Br), bond / ring / aromaticity deltas,
    molecular-weight deltas, ring-size topology, Morgan-FP
    Tanimoto similarity, functional-group deltas (OH, C=O, CH₃,
    amine), stereocenter count, and degree-of-unsaturation.

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
    """
    try:
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol is None:
            return np.full(NUM_DOMAIN_FEATURES, np.nan)

        reactant_mols = [Chem.MolFromSmiles(r)
                         for r in reactants_smiles.split(".")]

        # Product-side descriptors
        p_atoms = Counter(a.GetSymbol() for a in product_mol.GetAtoms())
        p_bonds = product_mol.GetNumBonds()
        p_rings = len(Chem.GetSSSR(product_mol))
        p_arom  = sum(1 for a in product_mol.GetAtoms() if a.GetIsAromatic())
        p_mw    = Descriptors.MolWt(product_mol)

        p_ring_sizes   = _ring_sizes(product_mol)
        p_heavy        = product_mol.GetNumHeavyAtoms()
        p_stereo       = len(Chem.FindMolChiralCenters(product_mol,
                                                        includeUnassigned=True))
        p_dbu = p_rings + sum(
            1 for b in product_mol.GetBonds()
            if b.GetBondTypeAsDouble() >= 2.0
        )

        # Reactant-side descriptors (aggregated across all valid reactant molecules)
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

        combined_reactant_mol = Chem.MolFromSmiles(reactants_smiles)

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

        if combined_reactant_mol is not None:
            tanimoto_sim = _tanimoto(product_mol, combined_reactant_mol)
        else:
            tanimoto_sim = 0.0

        heavy_atom_ratio = r_heavy / max(p_heavy, 1)

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

        stereocenter_delta = r_stereo - p_stereo
        unsaturation_delta = r_dbu - p_dbu

        return np.array([
            float(r_atoms.get('C', 0) - p_atoms.get('C', 0)),
            float(r_atoms.get('N', 0) - p_atoms.get('N', 0)),
            float(r_atoms.get('O', 0) - p_atoms.get('O', 0)),
            float(r_atoms.get('Cl', 0) - p_atoms.get('Cl', 0)),
            float(r_atoms.get('Br', 0) - p_atoms.get('Br', 0)),
            float(r_bonds - p_bonds),
            float(r_rings - p_rings),
            float(r_arom  - p_arom),
            float(r_mw    - p_mw),
            float(n_valid),
            float(sum(r_atoms.values()) - sum(p_atoms.values())),
            float(p_mw),
            float(r_mw),
            float(p_rings),
            float(r_rings),
            float(max_ring_delta),
            float(min_ring_delta),
            float(small_ring_count_r),
            float(ring_hist_dist),
            float(tanimoto_sim),
            float(heavy_atom_ratio),
            float(hydroxyl_delta),
            float(carbonyl_delta),
            float(methyl_delta),
            float(amine_delta),
            float(stereocenter_delta),
            float(unsaturation_delta),
        ])
    except Exception:
        return np.full(NUM_DOMAIN_FEATURES, np.nan)
