"""Domain feature extraction utilities for reaction-step featurization."""

import numpy as np
from collections import Counter
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import (
    CalcTPSA,
    CalcNumHBD,
    CalcNumHBA,
    CalcNumRotatableBonds,
)

NUM_DOMAIN_FEATURES = 39

# Internal Morgan fingerprint config used only for Tanimoto similarity
# inside the domain feature vector.  Kept independent of the user-facing
# fingerprint produced by ReactionStepFeaturizer so similarity is well
# defined regardless of which radius the featurizer is configured with.
_TANIMOTO_FP_RADIUS = 3
_TANIMOTO_FP_NBITS = 2048


def _morgan_fp(mol):
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=_TANIMOTO_FP_RADIUS, nBits=_TANIMOTO_FP_NBITS
    )


def extract_domain_features_single(
    product_smiles: str, reactants_smiles: str
) -> np.ndarray:
    """
    Extract hand-crafted domain features for one product-reactant pair.

    Returns 39 features.  Indices 0-14 (atom-count deltas C/N/O/Cl/Br,
    bond/ring/aromaticity/MW deltas, n_valid, total heavy delta, p_mw,
    r_mw, p_rings, r_rings) are preserved from the original 15-feature
    layout.  Indices 15-38 add chemistry-aware signals (S/F/P/I deltas,
    TPSA / HBD / HBA / rotbond / chiral deltas, max-/min-reactant MW,
    macrocycle flags, max-ring-size, Tanimoto similarities, heavy-atom
    ratio, ring-formed / ring-broken flags, plus absolute product TPSA /
    chiral / heavy-atom counts).

    Parameters
    ----------
    product_smiles : str
        SMILES of the target product.
    reactants_smiles : str
        SMILES of the proposed reactants (dot-separated when multiple).

    Returns
    -------
    features : np.ndarray, shape (NUM_DOMAIN_FEATURES,)
        1-D feature vector.  Returns a NaN vector on any parsing failure,
        so invalid rows are distinguishable from real data downstream.

    Examples
    --------
    >>> from deepretro.utils import extract_domain_features_single
    >>> feats = extract_domain_features_single("CCO", "CC.O")
    >>> feats.shape
    (39,)
    """
    try:
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol is None:
            return np.full(NUM_DOMAIN_FEATURES, np.nan)

        # Reactants may be multiple molecules joined by '.' (standard SMILES mixture notation)
        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]

        # Product-side descriptors
        p_atoms = Counter(
            a.GetSymbol() for a in product_mol.GetAtoms()
        )  # element frequency map
        p_bonds = product_mol.GetNumBonds()
        p_rings_obj = Chem.GetSSSR(product_mol)
        p_rings = len(p_rings_obj)
        p_arom = sum(1 for a in product_mol.GetAtoms() if a.GetIsAromatic())
        p_mw = Descriptors.MolWt(product_mol)
        p_tpsa = CalcTPSA(product_mol)
        p_hbd = CalcNumHBD(product_mol)
        p_hba = CalcNumHBA(product_mol)
        p_rotb = CalcNumRotatableBonds(product_mol)
        p_chiral = len(Chem.FindMolChiralCenters(product_mol, includeUnassigned=True))
        p_ring_sizes = [len(r) for r in p_rings_obj]
        p_max_ring = max(p_ring_sizes) if p_ring_sizes else 0
        p_is_macro = int(any(s > 8 for s in p_ring_sizes))
        p_n_heavy = product_mol.GetNumHeavyAtoms()

        # Reactant-side descriptors (aggregated across all valid reactant molecules)
        r_atoms, r_bonds, r_rings, r_arom, r_mw, n_valid = Counter(), 0, 0, 0, 0.0, 0
        r_tpsa = 0.0
        r_hbd = 0
        r_hba = 0
        r_rotb = 0
        r_chiral = 0
        r_ring_sizes_all: list[int] = []
        r_n_heavy = 0
        reactant_mws: list[float] = []
        for mol in reactant_mols:
            if mol:
                r_atoms += Counter(a.GetSymbol() for a in mol.GetAtoms())
                r_bonds += mol.GetNumBonds()
                mol_rings = Chem.GetSSSR(mol)
                r_rings += len(mol_rings)
                r_ring_sizes_all.extend(len(rg) for rg in mol_rings)
                r_arom += sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
                mol_mw = Descriptors.MolWt(mol)
                r_mw += mol_mw
                r_tpsa += CalcTPSA(mol)
                r_hbd += CalcNumHBD(mol)
                r_hba += CalcNumHBA(mol)
                r_rotb += CalcNumRotatableBonds(mol)
                r_chiral += len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
                r_n_heavy += mol.GetNumHeavyAtoms()
                reactant_mws.append(mol_mw)
                n_valid += 1  # count of successfully parsed reactant molecules

        r_max_ring = max(r_ring_sizes_all) if r_ring_sizes_all else 0
        r_is_macro = int(any(s > 8 for s in r_ring_sizes_all))

        # Tanimoto similarity (product vs first valid reactant; product vs best)
        valid_reactants = [m for m in reactant_mols if m is not None]
        if valid_reactants:
            p_fp = _morgan_fp(product_mol)
            tanimoto_first = DataStructs.TanimotoSimilarity(
                p_fp, _morgan_fp(valid_reactants[0])
            )
            tanimoto_max = max(
                DataStructs.TanimotoSimilarity(p_fp, _morgan_fp(m))
                for m in valid_reactants
            )
        else:
            tanimoto_first = 0.0
            tanimoto_max = 0.0

        max_reac_mw = max(reactant_mws) if reactant_mws else 0.0
        min_reac_mw = min(reactant_mws) if reactant_mws else 0.0

        # Heavy-atom ratio: real reactions usually have r_heavy >= p_heavy (some
        # heavy atoms come from leaving groups / reagents).
        heavy_ratio = (r_n_heavy / p_n_heavy) if p_n_heavy > 0 else 1.0
        ring_formed = int(p_rings > r_rings)
        ring_broken = int(p_rings < r_rings)

        return np.array(
            [
                # ── Indices 0-14: original 15 features (UNCHANGED order/values) ──
                float(r_atoms.get("C", 0) - p_atoms.get("C", 0)),  # 0  C delta
                float(r_atoms.get("N", 0) - p_atoms.get("N", 0)),  # 1  N delta
                float(r_atoms.get("O", 0) - p_atoms.get("O", 0)),  # 2  O delta
                float(r_atoms.get("Cl", 0) - p_atoms.get("Cl", 0)),  # 3  Cl delta
                float(r_atoms.get("Br", 0) - p_atoms.get("Br", 0)),  # 4  Br delta
                float(r_bonds - p_bonds),                            # 5  bond delta
                float(r_rings - p_rings),                            # 6  ring delta
                float(r_arom - p_arom),                              # 7  aromatic delta
                float(r_mw - p_mw),                                  # 8  mw delta
                float(n_valid),                                      # 9  n_valid
                float(sum(r_atoms.values()) - sum(p_atoms.values())),  # 10 heavy delta
                float(p_mw),                                         # 11 p_mw
                float(r_mw),                                         # 12 r_mw
                float(p_rings),                                      # 13 p_rings
                float(r_rings),                                      # 14 r_rings
                # ── Indices 15-38: new 24 features ──
                float(r_atoms.get("S", 0) - p_atoms.get("S", 0)),    # 15 S delta
                float(r_atoms.get("F", 0) - p_atoms.get("F", 0)),    # 16 F delta
                float(r_atoms.get("P", 0) - p_atoms.get("P", 0)),    # 17 P delta
                float(r_atoms.get("I", 0) - p_atoms.get("I", 0)),    # 18 I delta
                float(r_tpsa - p_tpsa),                              # 19 TPSA delta
                float(r_hbd - p_hbd),                                # 20 HBD delta
                float(r_hba - p_hba),                                # 21 HBA delta
                float(r_rotb - p_rotb),                              # 22 rotbond delta
                float(r_chiral - p_chiral),                          # 23 chiral delta
                float(max_reac_mw),                                  # 24 max reactant MW
                float(min_reac_mw),                                  # 25 min reactant MW
                float(p_is_macro),                                   # 26 p is_macrocycle
                float(r_is_macro),                                   # 27 r is_macrocycle
                float(p_max_ring),                                   # 28 p max ring size
                float(r_max_ring),                                   # 29 r max ring size
                float(p_max_ring - r_max_ring),                      # 30 delta max ring
                float(tanimoto_first),                               # 31 tanimoto first
                float(tanimoto_max),                                 # 32 tanimoto max
                float(heavy_ratio),                                  # 33 heavy ratio
                float(ring_formed),                                  # 34 ring_formed
                float(ring_broken),                                  # 35 ring_broken
                float(p_tpsa),                                       # 36 p_tpsa abs
                float(p_chiral),                                     # 37 p_n_chiral abs
                float(p_n_heavy),                                    # 38 p_n_atoms abs
            ]
        )
    except Exception:
        return np.full(NUM_DOMAIN_FEATURES, np.nan)
