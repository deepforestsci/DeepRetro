"""
Protecting group utilities for retrosynthesis workflows.

This module provides functions to mask protecting groups in SMILES strings for easier analysis.
"""
from rdkit import Chem
import re

# These are the protecting groups that we want to mask
PG_MAP = {
    "OMe": ("OC", "$"),
    "OBn": ("COCc1ccccc1", "%"),
    "OEt": ("COC", "&"),
}

def mask_protecting_groups_multisymbol(smiles: str) -> str:
    """
    Replace protecting groups in a SMILES string with single-character symbols.
    
    This function identifies and replaces common protecting groups in organic molecules
    with simplified single-character symbols to facilitate retrosynthetic analysis.
    The function first converts the input SMILES to canonical form, then replaces
    protecting groups in order of decreasing length to avoid partial matches.
    
    Parameters
    ----------
    smiles : str
        Input SMILES string representing an organic molecule. Can be any valid SMILES
        notation, including molecules with multiple disconnected components (e.g., salts).
        
    Returns
    -------
    str
        Modified SMILES string with protecting groups replaced by symbols:
        
        - 'OC' (OMe) → '$'
        - 'COCc1ccccc1' (OBn) → '%'
        - 'COC' (OEt) → '&'
        
        Returns 'INVALID_SMILES' if the input cannot be parsed by RDKit.
        Returns empty string if the input is an empty string.
    
    Examples
    --------
    >>> mask_protecting_groups_multisymbol("COC")  # doctest: +SKIP
    '&'
    >>> mask_protecting_groups_multisymbol("COCc1ccccc1")  # doctest: +SKIP
    '%'
    >>> mask_protecting_groups_multisymbol("CC(C)COCc1ccccc1")  # doctest: +SKIP
    'CC(C)%'
    >>> mask_protecting_groups_multisymbol("COC.COC")  # doctest: +SKIP
    '&'
    >>> mask_protecting_groups_multisymbol("invalid_smiles")  # doctest: +SKIP
    'INVALID_SMILES'
    >>> mask_protecting_groups_multisymbol("")  # doctest: +SKIP
    ''
    
    References
    ----------
    .. [1] Greene, T.W. and Wuts, P.G.M. (2006) Protective Groups in Organic Synthesis.
           4th Edition, John Wiley & Sons, Inc., Hoboken.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "INVALID_SMILES"

    full_smiles = Chem.MolToSmiles(mol, canonical=True)

    for pg_name, (pg_smiles, symbol) in sorted(PG_MAP.items(),
                                               key=lambda x: len(x[1][0]),
                                               reverse=True):
        full_smiles = full_smiles.replace(pg_smiles, symbol)

    # Clean up any edge-case SMILES formatting (e.g., ".%", "$.")
    full_smiles = re.sub(r"[\.\$%&]+", lambda m: m.group(0)[0], full_smiles)
    return full_smiles
