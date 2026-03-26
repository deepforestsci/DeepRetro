"""Domain feature extraction utilities for reaction-step featurization."""

import numpy as np
from collections import Counter
from rdkit import Chem
from rdkit.Chem import Descriptors

NUM_DOMAIN_FEATURES = 15


def extract_domain_features_single(
    product_smiles: str, reactants_smiles: str
) -> np.ndarray:
    """
    Extract hand-crafted domain features for one product-reactant pair.

    Computes atom-count deltas (C, N, O, Cl, Br), bond/ring/aromaticity
    deltas, molecular-weight deltas, and absolute counts.

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
    (15,)
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
        p_bonds = product_mol.GetNumBonds()  # total bond count
        p_rings = len(Chem.GetSSSR(product_mol))  # smallest set of smallest rings
        p_arom = sum(
            1 for a in product_mol.GetAtoms() if a.GetIsAromatic()
        )  # aromatic atom count
        p_mw = Descriptors.MolWt(product_mol)  # molecular weight (Da)

        # Reactant-side descriptors (aggregated across all valid reactant molecules)
        r_atoms, r_bonds, r_rings, r_arom, r_mw, n_valid = Counter(), 0, 0, 0, 0.0, 0
        for mol in reactant_mols:
            if mol:
                r_atoms += Counter(a.GetSymbol() for a in mol.GetAtoms())
                r_bonds += mol.GetNumBonds()
                r_rings += len(Chem.GetSSSR(mol))
                r_arom += sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
                r_mw += Descriptors.MolWt(mol)
                n_valid += 1  # count of successfully parsed reactant molecules

        # 15 features: deltas (reactant - product) capture what atoms/bonds are "consumed"
        return np.array(
            [
                float(r_atoms.get("C", 0) - p_atoms.get("C", 0)),  # carbon atom delta
                float(r_atoms.get("N", 0) - p_atoms.get("N", 0)),  # nitrogen atom delta
                float(r_atoms.get("O", 0) - p_atoms.get("O", 0)),  # oxygen atom delta
                float(
                    r_atoms.get("Cl", 0) - p_atoms.get("Cl", 0)
                ),  # chlorine atom delta
                float(
                    r_atoms.get("Br", 0) - p_atoms.get("Br", 0)
                ),  # bromine atom delta
                float(r_bonds - p_bonds),  # bond count delta
                float(r_rings - p_rings),  # ring count delta
                float(r_arom - p_arom),  # aromatic atom count delta
                float(r_mw - p_mw),  # molecular weight delta (Da)
                float(n_valid),  # number of valid reactant molecules
                float(
                    sum(r_atoms.values()) - sum(p_atoms.values())
                ),  # total heavy-atom delta
                float(p_mw),  # absolute product molecular weight
                float(r_mw),  # absolute reactant molecular weight (sum)
                float(p_rings),  # absolute product ring count
                float(r_rings),  # absolute reactant ring count (sum)
            ]
        )
    except Exception:
        return np.full(NUM_DOMAIN_FEATURES, np.nan)
