from rdkit import Chem
import re

# These are the protecting groups that we want to mask
PG_MAP = {
    "OMe": ("OC", "$"),
    "OBn": ("COCc1ccccc1", "%"),
    "OEt": ("COC", "&"),
}


def mask_protecting_groups_multisymbol(smiles: str) -> str:
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


if __name__ == "__main__":
    erythro_smiles = "CC(C)CC1=C(C(=O)OC)C(=C(C=C1OC)OC)C2OC(C(C(O2)COC(=O)C3=C(C(=C(C(=C3C(=O)OC)OC)OC)OC)OC)OC)OC"
    print(mask_protecting_groups_multisymbol(erythro_smiles))
