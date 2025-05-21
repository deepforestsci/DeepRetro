from src.utils.utils_molecule import is_valid_smiles
from src.utils.job_context import logger as context_logger
import os
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticCarbocycles, CalcNumAliphaticHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticRings, CalcNumAromaticCarbocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticHeterocycles, CalcNumAromaticRings
from rdkit.Chem.rdMolDescriptors import CalcNumBridgeheadAtoms

import rootutils

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

load_dotenv()

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
    log_message(
        f"Valid pathways: {valid_pathways} after stability check", logger)
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
                    if other_atom_idx in ring and bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bond_count += 0.5  # Count each end of double bond once

            # If number of π electrons is 4n (where n is positive integer)
            # Each double bond contributes 2 π electrons
            pi_electrons = int(double_bond_count * 2)
            if pi_electrons > 0 and pi_electrons % 4 == 0:
                results["issues"].append(
                    f"{len(ring)}-membered ring with {pi_electrons} π electrons (potential anti-aromatic)")

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
                        "Very large heterocycle with multiple heteroatoms (significant stability concern)"
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

    # ----------------- DETECT CARBOCATIONS -----------------
    # SMARTS patterns for detecting carbocations
    # sp2 carbocation pattern (trigonal planar)
    sp2_carbocation_pattern = Chem.MolFromSmarts("[C+;X3]")
    # sp carbocation pattern (linear)
    sp_carbocation_pattern = Chem.MolFromSmarts("[C+;X2]")
    # Non-stabilized carbocation (adjacent to electron-withdrawing groups)
    unstabilized_carbocation_pattern = Chem.MolFromSmarts(
        "[C+][F,Cl,Br,I,N+,S+,O+]")
    # Primary carbocation (less stable)
    primary_carbocation_pattern = Chem.MolFromSmarts("[CH2+][#6]")
    # Secondary carbocation (moderately stable)
    secondary_carbocation_pattern = Chem.MolFromSmarts("[CH+]([#6])[#6]")
    # Tertiary carbocation (more stable but still reactive)
    # tertiary_carbocation_pattern = Chem.MolFromSmarts("[C+]([#6])([#6])[#6]")
    # Allylic carbocation (stabilized by resonance)
    allylic_carbocation_pattern = Chem.MolFromSmarts("[C+]-[C]=[C]")
    # Benzylic carbocation (stabilized by resonance)
    benzylic_carbocation_pattern = Chem.MolFromSmarts("[C+]-c")

    # Check for each carbocation type
    if sp2_carbocation_pattern and mol.HasSubstructMatch(
            sp2_carbocation_pattern):
        matches = mol.GetSubstructMatches(sp2_carbocation_pattern)
        # Check if any of these are stabilized
        for match in matches:
            carbon_idx = match[0]
            carbon_atom = mol.GetAtomWithIdx(carbon_idx)
            neighbors = [
                mol.GetAtomWithIdx(a.GetIdx())
                for a in carbon_atom.GetNeighbors()
            ]

            # Check for stability factors
            is_allylic = mol.HasSubstructMatch(allylic_carbocation_pattern)
            is_benzylic = mol.HasSubstructMatch(benzylic_carbocation_pattern)

            if any(a.GetIsAromatic()
                   for a in neighbors) or is_allylic or is_benzylic:
                stability_score += 10
            else:
                results["issues"].append(
                    "Contains non-stabilized sp2 carbocation (highly unstable intermediate)"
                )
                stability_score -= 25  # Major penalty

    if sp_carbocation_pattern and mol.HasSubstructMatch(
            sp_carbocation_pattern):
        results["issues"].append(
            "Contains sp carbocation (highly unstable intermediate)")
        stability_score -= 30  # Major penalty

    if unstabilized_carbocation_pattern and mol.HasSubstructMatch(
            unstabilized_carbocation_pattern):
        results["issues"].append(
            "Contains carbocation adjacent to electron-withdrawing group (extremely unstable)"
        )
        stability_score -= 35  # Severe penalty

    # Check primary, secondary, tertiary carbocations
    if primary_carbocation_pattern and mol.HasSubstructMatch(
            primary_carbocation_pattern):
        results["issues"].append(
            "Contains primary carbocation (highly unstable)")
        stability_score -= 30
    elif secondary_carbocation_pattern and mol.HasSubstructMatch(
            secondary_carbocation_pattern):
        results["issues"].append(
            "Contains secondary carbocation (unstable intermediate)")
        stability_score -= 25

    # ----------------- DETECT CARBENES -----------------
    # Carbene patterns (neutral carbon with 2 bonds and no charge)
    singlet_carbene_pattern = Chem.MolFromSmarts("[C;X2;H0;+0]")
    # Carbenes with adjacent electron-withdrawing groups
    unstable_carbene_pattern = Chem.MolFromSmarts(
        "[C;X2;H0;+0][F,Cl,Br,I,N+,S+,O+]")
    # Carbenes in small rings (highly strained)
    ring_carbene_pattern_3 = Chem.MolFromSmarts("[C;X2;H0;+0]1[C,N,O][C,N,O]1")
    ring_carbene_pattern_4 = Chem.MolFromSmarts(
        "[C;X2;H0;+0]1[C,N,O][C,N,O][C,N,O]1")

    if singlet_carbene_pattern and mol.HasSubstructMatch(
            singlet_carbene_pattern):
        results["issues"].append(
            "Contains carbene (highly reactive intermediate)")
        stability_score -= 35  # Severe penalty

        # Check for additional destabilizing factors
        if unstable_carbene_pattern and mol.HasSubstructMatch(
                unstable_carbene_pattern):
            results["issues"].append(
                "Contains carbene adjacent to electron-withdrawing group (extremely unstable)"
            )
            stability_score -= 10  # Additional penalty

        if ring_carbene_pattern_3 and mol.HasSubstructMatch(
                ring_carbene_pattern_3):
            results["issues"].append(
                "Contains carbene in 3-membered ring (extremely unstable)")
            stability_score -= 15  # Additional penalty

        if ring_carbene_pattern_4 and mol.HasSubstructMatch(
                ring_carbene_pattern_4):
            results["issues"].append(
                "Contains carbene in 4-membered ring (highly unstable)")
            stability_score -= 10  # Additional penalty

    # ----------------- DETECT STRAINED 5-MEMBER RINGS WITH SMALL FUSED RINGS -----------------
    # Find 5-membered carbon rings with attached 3/4-membered rings containing heteroatoms
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if cyclopentane_pattern and mol.HasSubstructMatch(cyclopentane_pattern):
        # Find all 5-membered carbon rings
        five_mem_rings = [ring for ring in atom_rings if len(ring) == 5 and all(
            mol.GetAtomWithIdx(i).GetSymbol() == 'C' for i in ring)]

        # Find all 3/4-membered rings with heteroatoms
        small_hetero_rings = [ring for ring in atom_rings if len(ring) in [3, 4] and any(
            mol.GetAtomWithIdx(i).GetSymbol() != 'C' for i in ring)]

        # Check for fusion between these rings
        for five_ring in five_mem_rings:
            five_ring_set = set(five_ring)
            for small_ring in small_hetero_rings:
                small_ring_set = set(small_ring)
                shared_atoms = five_ring_set.intersection(small_ring_set)

                if len(shared_atoms) >= 1:
                    # Get heteroatom types in the small ring
                    heteroatoms = [mol.GetAtomWithIdx(i).GetSymbol()
                                   for i in small_ring if mol.GetAtomWithIdx(i).GetSymbol() != 'C']
                    unique_hetero = set(heteroatoms)

                    results["issues"].append(
                        f"5-membered carbon ring fused with {len(small_ring)}-membered ring containing {', '.join(unique_hetero)} (strained system)"
                    )
                    stability_score -= 40  # Significant penalty for this strained system

    # Penalize for extreme values of properties
    if abs(logp) > 10:
        stability_score -= 10
    if rotatable_bonds > 15:
        stability_score -= 10

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
