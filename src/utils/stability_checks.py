import os
import ast
from dotenv import load_dotenv
from rdkit import Chem
from litellm import completion
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticCarbocycles, CalcNumAliphaticHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticRings, CalcNumAromaticCarbocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticHeterocycles, CalcNumAromaticRings
from rdkit.Chem.rdMolDescriptors import CalcNumBridgeheadAtoms

import rootutils

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

load_dotenv()
from src.utils.job_context import logger as context_logger
from src.utils.utils_molecule import is_valid_smiles

ENABLE_LOGGING = False if os.getenv("ENABLE_LOGGING",
                                    "true").lower() == "false" else True


def log_message(message: str, logger=None):
    """Log the message

    Parameters
    ----------
    message : str
        The message to be logged
    logger : _type_, optional
        The logger object, by default None

    Returns
    -------
    None
    """
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def stability_checker(res_smiles: list):
    """Calls the LLM model to check the stability of a molecule

    Parameters
    ----------
    res_smiles : list
        list of SMILES strings of the molecule

    Returns
    -------
    int, list[str]
        Status code and list of stable SMILES
    """
    logger = context_logger.get() if ENABLE_LOGGING else None
    valid_pathways = []
    for idx, smile_list in enumerate(res_smiles):
        valid = []
        if isinstance(smile_list, list):
            for smiles in smile_list:
                if not is_valid_smiles(smiles):
                    log_message("Invalid SMILES string", logger)

                stability_dict = check_molecule_stability(smiles)
                log_message(f"Stability dict: {stability_dict}", logger)

                if stability_dict['assessment'] in [
                        "Likely stable", "Moderately stable"
                ]:
                    valid.append(smiles)
            if len(valid) == len(smile_list):
                valid_pathways.append(valid)
        else:
            if is_valid_smiles(smile_list):
                stability_dict = check_molecule_stability(smile_list)
                log_message(f"Stability dict: {stability_dict}", logger)

                if stability_dict['assessment'] in [
                        "Likely stable", "Moderately stable"
                ]:
                    valid_pathways.append([smile_list])
    log_message(f"Valid pathways: {valid_pathways} after stability check",
                logger)
    return 200, valid_pathways


def check_molecule_stability(smiles):
    """
    Performs heuristic checks on a molecule given its SMILES string to estimate stability.
    
    Parameters
    ----------
    smiles : str
        SMILES string of the molecule
        
    Returns
    -------
    results : dict
        Dictionary containing stability assessment and various metrics
    """
    # Initialize results dictionary
    results = {
        "valid_structure": False,
        "stability_score": 0,
        "issues": [],
        "metrics": {},
        "ring_data": {},
        "atom_data": {},
        "assessment": "",
    }

    # Parse SMILES and check if it's a valid structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        results["issues"].append("Invalid SMILES string or cannot be parsed")
        return results

    results["valid_structure"] = True

    # Calculate basic properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    results["metrics"] = {
        "molecular_weight": mw,
        "logP": logp,
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": rotatable_bonds
    }

    # Get detailed ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()

    # Get atom and bond counts
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    num_aromatic_atoms = len(mol.GetAromaticAtoms())
    num_aliphatic_atoms = num_heavy_atoms - num_aromatic_atoms

    # Calculate ring-related descriptors
    num_aliphatic_carbocycles = CalcNumAliphaticCarbocycles(mol)
    num_aliphatic_heterocycles = CalcNumAliphaticHeterocycles(mol)
    num_aliphatic_rings = CalcNumAliphaticRings(mol)
    num_aromatic_carbocycles = CalcNumAromaticCarbocycles(mol)
    num_aromatic_heterocycles = CalcNumAromaticHeterocycles(mol)
    num_aromatic_rings = CalcNumAromaticRings(mol)
    num_bridgehead_atoms = CalcNumBridgeheadAtoms(mol)

    # Store ring data
    results["ring_data"] = {
        "num_rings": len(atom_rings),
        "atom_rings": [list(ring) for ring in atom_rings],
        "bond_rings": [list(ring) for ring in bond_rings],
        "num_aliphatic_carbocycles": num_aliphatic_carbocycles,
        "num_aliphatic_heterocycles": num_aliphatic_heterocycles,
        "num_aliphatic_rings": num_aliphatic_rings,
        "num_aromatic_carbocycles": num_aromatic_carbocycles,
        "num_aromatic_heterocycles": num_aromatic_heterocycles,
        "num_aromatic_rings": num_aromatic_rings,
        "num_bridgehead_atoms": num_bridgehead_atoms
    }

    # Store atom data
    results["atom_data"] = {
        "num_atoms": num_atoms,
        "num_bonds": num_bonds,
        "num_heavy_atoms": num_heavy_atoms,
        "num_heavy_bonds":
        num_bonds,  # Same as num_bonds since RDKit doesn't count hydrogen bonds separately
        "num_aromatic_atoms": num_aromatic_atoms,
        "num_aliphatic_atoms": num_aliphatic_atoms
    }

    # Check for strained rings (small rings are often unstable)
    for ring in atom_rings:
        if len(ring) < 3:
            results["issues"].append(
                f"Highly strained ring of size {len(ring)}")
        elif len(ring) == 3 and any(
                mol.GetAtomWithIdx(i).GetSymbol() != 'C' for i in ring):
            results["issues"].append(
                "Three-membered heterocycle (potentially unstable)")
        elif len(ring) == 4 and any(
                mol.GetAtomWithIdx(i).GetSymbol() != 'C' for i in ring):
            results["issues"].append(
                "Four-membered heterocycle (potentially unstable)")

    # ----------------- DETECT ANTI-AROMATIC COMPOUNDS -----------------
    # Check for anti-aromatic systems (4n π electrons)
    # Common anti-aromatic patterns include cyclobutadiene, cyclooctatetraene, etc.
    patt_cyclobutadiene = Chem.MolFromSmarts(
        "c1ccc1")  # 4-membered fully conjugated ring
    patt_cyclooctatetraene = Chem.MolFromSmarts(
        "C1=CC=CC=CC=C1")  # Non-planar COT
    patt_pentalene = Chem.MolFromSmarts("c1cc2cccc2c1")  # Pentalene pattern

    anti_aromatic_patterns = [
        (patt_cyclobutadiene, "cyclobutadiene-like (anti-aromatic)"),
        (patt_cyclooctatetraene,
         "cyclooctatetraene-like (potential anti-aromatic)"),
        (patt_pentalene, "pentalene-like (potential anti-aromatic)")
    ]

    for patt, name in anti_aromatic_patterns:
        if patt and mol.HasSubstructMatch(patt):
            results["issues"].append(f"Contains {name} motif")

    # Also detect rings with 4n π electrons
    for ring in atom_rings:
        # Only consider rings of sizes that could be anti-aromatic
        if len(ring) in [4, 8, 12, 16]:
            # Count conjugated double bonds in the ring
            double_bond_count = 0
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]

            for atom in ring_atoms:
                if atom.GetIsAromatic():
                    double_bond_count += 0.5

                for bond in atom.GetBonds():
                    # Check if bond is in this ring and is a double bond
                    other_atom_idx = bond.GetOtherAtomIdx(atom.GetIdx())
                    if other_atom_idx in ring and bond.GetBondType(
                    ) == Chem.BondType.DOUBLE:
                        double_bond_count += 0.5  # Count each end of double bond once

            # If number of π electrons is 4n (where n is positive integer)
            pi_electrons = int(double_bond_count *
                               2)  # Each double bond contributes 2 π electrons
            if pi_electrons > 0 and pi_electrons % 4 == 0:
                results["issues"].append(
                    f"{len(ring)}-membered ring with {pi_electrons} π electrons (potential anti-aromatic)"
                )

    # ----------------- DETECT FUSED 3-4 MEMBERED RINGS -----------------
    # Look for fused small rings (which create highly strained systems)
    small_rings = []
    for ring in atom_rings:
        if len(ring) <= 4:
            small_rings.append(set(ring))

    # Check for fused small rings (rings that share atoms)
    fused_small_rings_detected = False
    for i in range(len(small_rings)):
        for j in range(i + 1, len(small_rings)):
            shared_atoms = small_rings[i].intersection(small_rings[j])
            if len(shared_atoms) > 0:
                fused_small_rings_detected = True
                ring1_size = len(small_rings[i])
                ring2_size = len(small_rings[j])
                results["issues"].append(
                    f"Fused {ring1_size} and {ring2_size}-membered rings (highly strained system)"
                )

    if fused_small_rings_detected:
        # This is a serious stability concern, add an explicit warning
        results["issues"].append(
            "WARNING: Fused small rings create highly strained and potentially explosive compounds"
        )

    # ----------------- DETECT LARGE HETEROCYCLES -----------------
    # Large heterocycles (>7 members) can be unstable
    for ring in atom_rings:
        if len(ring) >= 7:
            # Check if ring contains heteroatoms
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            heteroatoms = [
                atom for atom in ring_atoms
                if atom.GetSymbol() not in ['C', 'H']
            ]

            if len(heteroatoms) > 0:
                hetero_symbols = [atom.GetSymbol() for atom in heteroatoms]
                unique_hetero = set(hetero_symbols)
                results["issues"].append(
                    f"Large ({len(ring)}-membered) heterocycle with {', '.join(unique_hetero)} (potentially unstable)"
                )

                # More severe warning for very large heterocycles (>10 members)
                if len(ring) > 10 and len(heteroatoms) >= 3:
                    results["issues"].append(
                        f"Very large heterocycle with multiple heteroatoms (significant stability concern)"
                    )

    # Assess complex ring systems
    if num_bridgehead_atoms > 0:
        # Complex polycyclic systems can be strained
        if any(len(ring) <= 4
               for ring in atom_rings) and num_bridgehead_atoms >= 2:
            results["issues"].append(
                "Strained polycyclic system with bridgehead atoms")

    # Calculate stability score (0-100, higher is more stable)
    # Start with 100 and subtract for various issues
    stability_score = 100

    # Penalize for each identified issue
    stability_score -= len(results["issues"]) * 15

    # Penalize for extreme values of properties
    if abs(logp) > 10: stability_score -= 10
    if rotatable_bonds > 15: stability_score -= 10

    # Assess stability based on ring structure
    if num_aromatic_rings > 0:
        # Aromatic rings typically enhance stability
        stability_score += min(num_aromatic_rings * 5, 15)  # Cap bonus at 15

    if num_aliphatic_heterocycles > 0 and num_aliphatic_rings <= 3:
        # Small number of aliphatic heterocycles can be unstable
        stability_score -= 5 * num_aliphatic_heterocycles

    # Penalize for complex strained systems
    if num_bridgehead_atoms > 0 and any(len(ring) <= 5 for ring in atom_rings):
        stability_score -= num_bridgehead_atoms * 5

    # Additional penalties for new detected issues
    if fused_small_rings_detected:
        stability_score -= 30  # Severe penalty for fused small rings

    # Penalize anti-aromatic structures
    for issue in results["issues"]:
        if "anti-aromatic" in issue:
            stability_score -= 25
        elif "π electrons" in issue:
            stability_score -= 20

    # Penalize large heterocycles
    large_heterocycle_count = sum(
        1 for issue in results["issues"]
        if "heterocycle" in issue and "large" in issue.lower())
    if large_heterocycle_count > 0:
        stability_score -= large_heterocycle_count * 10

    # Cap the score between 0 and 100
    stability_score = max(0, min(100, stability_score))
    results["stability_score"] = stability_score

    # Overall stability assessment
    if stability_score >= 80:
        results["assessment"] = "Likely stable"
    elif stability_score >= 50:
        results["assessment"] = "Moderately stable"
    else:
        results["assessment"] = "Potentially unstable"

    return results
