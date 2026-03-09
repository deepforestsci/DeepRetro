import os

import rootutils
import structlog
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

logger = structlog.get_logger()

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

RXN_CLASSIFICATION_MODEL_PATH = (
    f"{root_dir}/{os.getenv('RXN_CLASSIFICATION_MODEL_PATH')}")


def is_valid_smiles(smiles: str) -> bool:
    """Check if the SMILES string is valid

    Parameters
    ----------
    smiles : str
        smiles string

    Returns
    -------
    bool
        True if the smiles is valid, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return False
    if mol is None:
        return False
    return True


def substructure_matching(target_smiles: str, query_smiles: str) -> int:
    """Check if the query substructure is present in the target molecule

    Parameters
    ----------
    target_smiles : str
        SMILES string of the target molecule
    query_smiles : str
        SMILES string of the query molecule

    Returns
    -------
    int
        1 if the query substructure is present in the target molecule, 0 otherwise
    """
    try:
        target_molecule = Chem.MolFromSmiles(target_smiles)
    except Exception:
        logger.error("Error parsing target molecule", smiles=target_smiles)

    try:
        query_molecule = Chem.MolFromSmiles(query_smiles)
    except Exception:
        logger.error("Error parsing query molecule", smiles=query_smiles)

    # Check if the query substructure is present in the target molecule
    try:
        if target_molecule.HasSubstructMatch(query_molecule):
            return 1
        else:
            return 0
    except Exception:
        return 0


def validity_check(molecule, res_molecules, res_explanations, res_confidence):
    """Check the validity of the molecules obtained from LLM

    Parameters
    ----------
    molecule : str
        Target molecule for retrosynthesis
    res_molecules : list
        List of molecules obtained from LLM
    res_explanations : list
        List of explanations obtained from LLM
    res_confidence : list
        List of confidence scores obtained from LLM

    Returns
    -------
    list
        List of valid pathways
    list
        List of valid explanations
    list
        List of valid confidence scores
    """
    valid_pathways = []
    valid_explanations = []
    valid_confidence = []
    for idx, smile_list in enumerate(res_molecules):
        valid = []
        if isinstance(smile_list, list):
            for smiles in smile_list:
                if is_valid_smiles(smiles):
                    if are_molecules_same(molecule, smiles):
                        logger.warning("Molecule is same as target",
                                       molecule=molecule,
                                       smiles=smiles)
                    elif substructure_matching(smiles, molecule):
                        logger.warning("Molecule is substructure of target",
                                       molecule=molecule,
                                       smiles=smiles)
                    else:
                        valid.append(smiles)
                else:
                    logger.warning("Invalid SMILES in pathway",
                                   molecule=molecule,
                                   smiles=smiles)
            if len(valid) == len(smile_list):
                valid_pathways.append(valid)
                valid_explanations.append(res_explanations[idx])
                valid_confidence.append(res_confidence[idx])
        else:
            if is_valid_smiles(smile_list):
                if are_molecules_same(molecule, smile_list):
                    logger.warning("Molecule is same as target",
                                   molecule=molecule,
                                   smiles=smile_list)
                elif substructure_matching(smile_list, molecule):
                    logger.warning("Molecule is substructure of target",
                                   molecule=molecule,
                                   smiles=smile_list)
                else:
                    valid_pathways.append([smile_list])
                    valid_explanations.append(res_explanations[idx])
                    valid_confidence.append(res_confidence[idx])
            else:
                logger.warning("Invalid SMILES", smiles=smile_list)
    logger.info("Validity check complete", valid_pathways=len(valid_pathways))
    return valid_pathways, valid_explanations, valid_confidence


def calc_mol_wt(mol: str) -> float:
    """Calculate the molecular weight of a molecule

    Parameters
    ----------
    mol : str
        SMILES string of the molecule

    Returns
    -------
    float
        molecular weight of the molecule
    """
    try:
        mol_wt = ExactMolWt(Chem.MolFromSmiles(mol))
    except Exception:
        mol_wt = 0.0
        logger.error("Error calculating molecular weight", smiles=mol)
    return mol_wt


def calc_chemical_formula(mol: str):
    """Calculate the chemical formula of a molecule

    Parameters
    ----------
    mol : str
        SMILES string of the molecule

    Returns
    -------
    str
        molecular formula of the molecule
    """
    try:
        formula = CalcMolFormula(Chem.MolFromSmiles(mol))
    except Exception:
        formula = "N/A"
        logger.error("Error calculating chemical formula", smiles=mol)
    return formula


def are_molecules_same(smiles1: str, smiles2: str) -> bool:
    """Check if two molecules are the same

    Parameters
    ----------
    smiles1 : str
        SMILES string of the first molecule
    smiles2 : str
        SMILES string of the second molecule

    Returns
    -------
    bool
        True if the molecules are the same, False otherwise
    """
    # Convert SMILES strings to RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        raise ValueError("Invalid SMILES string provided.")

    # Get canonical SMILES for both molecules
    canonical_smiles1 = Chem.MolToSmiles(mol1, canonical=True)
    canonical_smiles2 = Chem.MolToSmiles(mol2, canonical=True)

    # Alternatively, compare molecular fingerprints
    fingerprint1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1,
                                                                  radius=2,
                                                                  nBits=1024)
    fingerprint2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2,
                                                                  radius=2,
                                                                  nBits=1024)

    # Check if canonical SMILES or fingerprints match
    if canonical_smiles1 == canonical_smiles2:
        return True
    elif fingerprint1 == fingerprint2:
        return True
    else:
        return False


def compute_fingerprint(smiles, radius=2, nBits=2048):
    """Compute the fingerprint of a molecule

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule
    radius : int, optional
    nBits : int, optional
        Number of bits in the fingerprint

    Returns
    -------
    list
        Fingerprint of the molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,
                                                        radius,
                                                        nBits=nBits)
    return list(fingerprint)


def detect_seven_member_rings(smiles) -> bool:
    """
    Detects 7-member rings in a molecule given its SMILES string.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    bool
        True if 7-member rings are present, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")

    # Retrieve ring information as tuples of atom indices.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Filter rings by the number of atoms.
    rings_7 = [ring for ring in atom_rings if len(ring) == 7]

    if len(rings_7) > 0:
        return True
    return False


def detect_eight_member_rings(smiles) -> bool:
    """
    Detects 8-member rings in a molecule given its SMILES string.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    bool
        True if 8-member rings are present, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")

    # Retrieve ring information as tuples of atom indices.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Filter rings by the number of atoms.
    rings_8 = [ring for ring in atom_rings if len(ring) == 8]

    if len(rings_8) > 0:
        return True
    return False
