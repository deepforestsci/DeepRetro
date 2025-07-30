from rdkit import Chem
import re

# These are the protecting groups that we want to unmask (reverse of PG_MAP)
DEPROTECT_MAP = {
    "$": "OC",  # OMe
    "%": "COCc1ccccc1",  # OBn
    "&": "COC",  # OEt
}


def unmask_protecting_groups_multisymbol(smiles: str) -> str:
    """
    Replace single-character symbols in a SMILES string with their original protecting groups.
    
    This function reverses the masking process by converting simplified single-character
    symbols back to their original protecting group SMILES representations. This is
    useful for converting masked molecules back to their full chemical representations
    after retrosynthetic analysis.
    
    Parameters
    ----------
    smiles : str
        Input SMILES string containing masked protecting group symbols. Can be any
        valid SMILES notation with symbols '$', '%', '&' representing protecting groups.
        
    Returns
    -------
    str
        Modified SMILES string with symbols replaced by original protecting groups:
        
        - '$' → 'OC' (OMe)
        - '%' → 'COCc1ccccc1' (OBn)  
        - '&' → 'COC' (OEt)
        
        Returns 'INVALID_SMILES' if the input cannot be parsed by RDKit.
        Returns empty string if the input is an empty string.
    
    Examples
    --------
    >>> unmask_protecting_groups_multisymbol("&")
    'COC'
    
    >>> unmask_protecting_groups_multisymbol("%")
    'COCc1ccccc1'
    
    >>> unmask_protecting_groups_multisymbol("CC(C)%")
    'CC(C)COCc1ccccc1'
    
    >>> unmask_protecting_groups_multisymbol("&.&")
    'COC.COC'
    
    >>> unmask_protecting_groups_multisymbol("CC(C)$%&")
    'CC(C)OCCOCc1ccccc1COC'
    
    >>> unmask_protecting_groups_multisymbol("invalid_smiles")
    'INVALID_SMILES'
    
    >>> unmask_protecting_groups_multisymbol("")
    ''
      
    References
    ----------
    .. [1] Greene, T.W. and Wuts, P.G.M. (2006) Protective Groups in Organic Synthesis.
           4th Edition, John Wiley & Sons, Inc., Hoboken.
    """
    if not smiles:
        return ""

    # Handle special case for "INVALID_SMILES"
    if smiles == "INVALID_SMILES":
        return "INVALID_SMILES"

    # First, replace all symbols with their corresponding protecting groups
    unmasked_smiles = smiles
    for symbol, pg_smiles in DEPROTECT_MAP.items():
        unmasked_smiles = unmasked_smiles.replace(symbol, pg_smiles)

    # Validate that the resulting SMILES is valid
    mol = Chem.MolFromSmiles(unmasked_smiles)
    if mol is None:
        return "INVALID_SMILES"

    # Return canonical form of the unmasked SMILES
    return Chem.MolToSmiles(mol, canonical=True)


def get_protecting_group_info() -> dict:
    """
    Get information about the protecting groups used in the masking/unmasking process.
    
    Returns
    -------
    dict
        Dictionary containing information about each protecting group:
        - symbol: The single-character symbol used for masking
        - smiles: The SMILES representation of the protecting group
        - name: The common name of the protecting group
        - description: Brief description of the protecting group
    """
    return {
        "$": {
            "smiles": "OC",
            "name": "OMe",
            "description": "Methoxy protecting group"
        },
        "%": {
            "smiles": "COCc1ccccc1",
            "name": "OBn",
            "description": "Benzyl protecting group"
        },
        "&": {
            "smiles": "COC",
            "name": "OEt",
            "description": "Ethoxy protecting group"
        }
    }


def validate_masked_smiles(smiles: str) -> bool:
    """
    Validate if a SMILES string contains valid masking symbols.
    
    Parameters
    ----------
    smiles : str
        Input SMILES string to validate
        
    Returns
    -------
    bool
        True if the SMILES contains valid masking symbols ('$', '%', '&')
        and other valid SMILES characters, False otherwise
    """
    if not smiles:
        return True

    # Check if the SMILES contains any of our masking symbols
    has_masking_symbols = any(symbol in smiles
                              for symbol in DEPROTECT_MAP.keys())

    if not has_masking_symbols:
        return False

    # Define valid characters for masked SMILES (including our symbols)
    valid_chars = set("$%&") | set(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]{}@+-=#$%&*./\\"
    )

    # Check if all characters in the SMILES are valid
    return all(char in valid_chars for char in smiles)
