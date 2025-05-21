import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt

import joblib
import rootutils
from src.variables import REACTION_ENCODING_NAMES, ENCODING_SCALABILITY
from src.cache import cache_results
from src.utils.job_context import logger as context_logger

root_dir = rootutils.setup_root(
    __file__, indicator=".project-root", pythonpath=True)

RXN_CLASSIFICATION_MODEL_PATH = f"{root_dir}/{os.getenv('RXN_CLASSIFICATION_MODEL_PATH')}"
ENABLE_LOGGING = False if os.getenv(
    "ENABLE_LOGGING", "true").lower() == "false" else True


def log_message(message: str, logger=None):
    """Log the message"""
    if logger is not None:
        log_message(message)
    else:
        print(message)


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
    logger = context_logger.get() if ENABLE_LOGGING else None

    # Convert SMILES to RDKit molecule objects
    try:
        target_molecule = Chem.MolFromSmiles(target_smiles)
    except Exception:
        log_message(
            f"Error in parsing target molecule: {target_smiles}", logger)
        return 0

    try:
        query_molecule = Chem.MolFromSmiles(query_smiles)
    except Exception:
        log_message(f"Error in parsing query molecule: {query_smiles}", logger)
        return 0

    # Check if the query substructure is present in the target molecule
    try:
        if target_molecule.HasSubstructMatch(query_molecule):
            return 1
        else:
            return 0
    except Exception:
        return 0


@cache_results
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
    logger = context_logger.get() if ENABLE_LOGGING else None
    valid_pathways = []
    valid_explanations = []
    valid_confidence = []
    for idx, smile_list in enumerate(res_molecules):
        valid = []
        if isinstance(smile_list, list):
            for smiles in smile_list:
                if is_valid_smiles(smiles):
                    if are_molecules_same(molecule, smiles):
                        log_message(
                            f"Molecule : {molecule} is same as target molecule",
                            logger)
                    elif substructure_matching(smiles, molecule):
                        log_message(
                            f"Molecule : {molecule} is substructure of target molecule",
                            logger)
                    else:
                        valid.append(smiles)
                else:
                    log_message(
                        f"Molecule : {molecule} is invalid or cannot be parsed",
                        logger)
            if len(valid) == len(smile_list):
                valid_pathways.append(valid)
                valid_explanations.append(res_explanations[idx])
                valid_confidence.append(res_confidence[idx])
        else:
            if is_valid_smiles(smile_list):
                if are_molecules_same(molecule, smile_list):
                    log_message("Molecule is same as target molecule", logger)
                elif substructure_matching(smile_list, molecule):
                    log_message(
                        f"Molecule : {molecule} is substructure of target molecule {smile_list}",
                        logger)
                else:
                    valid_pathways.append([smile_list])
                    valid_explanations.append(res_explanations[idx])
                    valid_confidence.append(res_confidence[idx])
            else:
                log_message("Molecule is invalid or cannot be parsed", logger)
    log_message(
        f"Obtained {len(valid_pathways)} valid pathways after validity test: {valid_pathways}",
        logger)
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
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        mol_wt = ExactMolWt(Chem.MolFromSmiles(mol))
    except Exception:
        mol_wt = 0.0
        log_message(f"Error in calculating molecular weight: {mol}", logger)
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
    logger = context_logger.get() if ENABLE_LOGGING else None
    try:
        formula = CalcMolFormula(Chem.MolFromSmiles(mol))
    except Exception:
        formula = "N/A"
        log_message(f"Error in calculating formula: {mol}", logger)
    return formula


def are_molecules_same(smiles1: str, smiles2: str) -> bool:
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,
                                                        radius,
                                                        nBits=nBits)
    return list(fingerprint)


def sub_structure_matching(target_smiles: str, query_smiles: str) -> bool:
    """Check if the query substructure is present in the target molecule"""
    target_molecule = Chem.MolFromSmiles(target_smiles)
    query_molecule = Chem.MolFromSmiles(query_smiles)

    if target_molecule.HasSubstructMatch(query_molecule):
        return True
    else:
        return False


def get_reaction_type(mol1, mol2, model_path):
    """Get the reaction type of a reaction"""
    clf = joblib.load(model_path)
    mol1_fingerprint = compute_fingerprint(mol1)
    mol2_fingerprint = compute_fingerprint(mol2)
    reaction_type = clf.predict([mol1_fingerprint + mol2_fingerprint])
    return REACTION_ENCODING_NAMES[reaction_type[0]], reaction_type[0]


def calc_confidence_estimate(probability: float) -> float:
    """Calculate the confidence estimate based on the probability

    Parameters
    ----------
    probability : float
        Probability of the prediction

    Returns
    -------
    float
        Confidence estimate
    """
    if isinstance(probability, list):
        probability = probability[0]
    if probability < 0.3:
        probability = 1 - probability
    elif probability < 0.45 and probability >= 0.3:
        probability += 0.5
    elif probability < 0.6 and probability >= 0.45:
        probability += 0.3

    # limit the confidence estimate to 2 decimal places, round to the
    # nearest 0.01
    probability = round(probability, 2)
    if probability > 0.99:
        probability = 0.99
    return probability


def calc_scalability_index(mol1, mol2):
    """Calculate the scalability index of a reaction"""
    _, type = get_reaction_type(mol1, mol2, RXN_CLASSIFICATION_MODEL_PATH)
    return str(ENCODING_SCALABILITY[type])


def calc_yield(mol1, mol2):
    """Calculate the yield of a reaction"""
    return "#"


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
