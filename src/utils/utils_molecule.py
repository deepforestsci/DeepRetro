from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
import logging
import joblib
from src.variables import REACTION_ENCODING_NAMES
from src.cache import cache_results

logger = logging.getLogger(__name__)


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
    except:
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
    # Convert SMILES to RDKit molecule objects
    try:
        target_molecule = Chem.MolFromSmiles(target_smiles)
    except:
        logger.info(f"Error in parsing target molecule: {target_smiles}")

    try:
        query_molecule = Chem.MolFromSmiles(query_smiles)
    except:
        logger.info(f"Error in parsing query molecule: {query_smiles}")

    # Check if the query substructure is present in the target molecule
    try:
        if target_molecule.HasSubstructMatch(query_molecule):
            return 1
        else:
            return 0
    except:
        return 0


@cache_results
def validity_check(molecule, res_molecules, res_explanations, res_confidence):
    valid_pathways = []
    valid_explanations = []
    valid_confidence = []
    for idx, smile_list in enumerate(res_molecules):
        valid = []
        if isinstance(smile_list, list):
            for smiles in smile_list:
                if is_valid_smiles(smiles):
                    if are_molecules_same(molecule, smiles):
                        logger.info(
                            f"Molecule : {molecule} is same as target molecule"
                        )
                    elif substructure_matching(smiles, molecule):
                        logger.info(
                            f"Molecule : {molecule} is substructure of target molecule"
                        )
                    else:
                        valid.append(smiles)
                else:
                    logger.info(
                        f"Molecule : {molecule} is invalid or cannot be parsed"
                    )
            if len(valid) >= 2:
                valid_pathways.append(valid)
                valid_explanations.append(res_explanations[idx])
                valid_confidence.append(res_confidence[idx])
        else:
            if is_valid_smiles(smile_list):
                if are_molecules_same(molecule, smiles):
                    logger.info("Molecule is same as target molecule")
                elif substructure_matching(smiles, molecule):
                    logger.info(
                        f"Molecule : {molecule} is substructure of target molecule {smiles}"
                    )
                else:
                    valid_pathways.append([smile_list])
                    valid_explanations.append(res_explanations[idx])
                    valid_confidence.append(res_confidence[idx])
            else:
                logger.info("Molecule is invalid or cannot be parsed")
    logger.info(
        f"Obtained {len(valid_pathways)} valid pathways after validity test: {valid_pathways}"
    )
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
    except:
        mol_wt = 0.0
        logger.info(f"Error in calculating molecular weight: {mol}")
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
    except:
        formula = "N/A"
        logger.info(f"Error in calculating formula: {mol}")
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
