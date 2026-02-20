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
        1-D feature vector.  Returns zeros on any parsing failure.
    """
    try:
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol is None:
            return np.zeros(NUM_DOMAIN_FEATURES)

        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]

        p_atoms = Counter(a.GetSymbol() for a in product_mol.GetAtoms())
        p_bonds = product_mol.GetNumBonds()
        p_rings = len(Chem.GetSSSR(product_mol))
        p_arom = sum(1 for a in product_mol.GetAtoms() if a.GetIsAromatic())
        p_mw = Descriptors.MolWt(product_mol)  # type: ignore

        r_atoms, r_bonds, r_rings, r_arom, r_mw, n_valid = Counter(), 0, 0, 0, 0.0, 0
        for mol in reactant_mols:
            if mol:
                r_atoms += Counter(a.GetSymbol() for a in mol.GetAtoms())
                r_bonds += mol.GetNumBonds()
                r_rings += len(Chem.GetSSSR(mol))
                r_arom += sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
                r_mw += Descriptors.MolWt(mol)  # type: ignore
                n_valid += 1

        return np.array(
            [
                float(r_atoms.get("C", 0) - p_atoms.get("C", 0)),
                float(r_atoms.get("N", 0) - p_atoms.get("N", 0)),
                float(r_atoms.get("O", 0) - p_atoms.get("O", 0)),
                float(r_atoms.get("Cl", 0) - p_atoms.get("Cl", 0)),
                float(r_atoms.get("Br", 0) - p_atoms.get("Br", 0)),
                float(r_bonds - p_bonds),
                float(r_rings - p_rings),
                float(r_arom - p_arom),
                float(r_mw - p_mw),
                float(n_valid),
                float(sum(r_atoms.values()) - sum(p_atoms.values())),
                float(p_mw),
                float(r_mw),
                float(p_rings),
                float(r_rings),
            ]
        )
    except Exception:
        return np.zeros(NUM_DOMAIN_FEATURES)
