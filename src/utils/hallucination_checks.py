import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops
from collections import Counter
import numpy as np


def compare_molecules(reactant_smiles, product_smiles):
    """
    Compare reactant and product molecules to detect potential hallucinations or inconsistencies.
    
    Parameters:
    -----------
    reactant_smiles : str
        SMILES string of the reactant molecule
    product_smiles : str
        SMILES string of the product molecule
        
    Returns:
    --------
    dict
        Dictionary containing validation results and detected issues
    """
    results = {
        "valid_reactant": False,
        "valid_product": False,
        "atom_count_consistent": False,
        "ring_size_changes": [],
        "substituent_position_changes": [],
        "detected_issues": []
    }

    # Check if SMILES strings are valid
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    if reactant_mol is None:
        results["detected_issues"].append("Invalid reactant SMILES string")
        return results
    else:
        results["valid_reactant"] = True

    if product_mol is None:
        results["detected_issues"].append("Invalid product SMILES string")
        return results
    else:
        results["valid_product"] = True

    # Get basic molecule properties
    reactant_atoms = Counter(
        [atom.GetSymbol() for atom in reactant_mol.GetAtoms()])
    product_atoms = Counter(
        [atom.GetSymbol() for atom in product_mol.GetAtoms()])

    # Check atom count consistency
    for atom_symbol in set(
            list(reactant_atoms.keys()) + list(product_atoms.keys())):
        if reactant_atoms.get(atom_symbol,
                              0) != product_atoms.get(atom_symbol, 0):
            results["detected_issues"].append(
                f"Atom count mismatch for {atom_symbol}: "
                f"Reactant has {reactant_atoms.get(atom_symbol, 0)}, "
                f"Product has {product_atoms.get(atom_symbol, 0)}")

    if not any("Atom count mismatch" in issue
               for issue in results["detected_issues"]):
        results["atom_count_consistent"] = True

    # Check for ring size changes
    reactant_rings = Chem.GetSSSR(reactant_mol)
    product_rings = Chem.GetSSSR(product_mol)

    reactant_ring_sizes = [len(ring) for ring in reactant_rings]
    product_ring_sizes = [len(ring) for ring in product_rings]

    # Sort ring sizes for easier comparison
    reactant_ring_sizes.sort()
    product_ring_sizes.sort()

    if reactant_ring_sizes != product_ring_sizes:
        results["detected_issues"].append(
            f"Ring size change detected: Reactant rings {reactant_ring_sizes}, "
            f"Product rings {product_ring_sizes}")

        # Report specific ring changes
        for r_size in reactant_ring_sizes:
            if reactant_ring_sizes.count(r_size) > product_ring_sizes.count(
                    r_size):
                results["ring_size_changes"].append(
                    f"{r_size}-membered ring removed")

        for p_size in product_ring_sizes:
            if product_ring_sizes.count(p_size) > reactant_ring_sizes.count(
                    p_size):
                results["ring_size_changes"].append(
                    f"{p_size}-membered ring added")

    # Check for aromatic ring changes
    reactant_aromatic_atoms = set([
        atom.GetIdx() for atom in reactant_mol.GetAtoms()
        if atom.GetIsAromatic()
    ])
    product_aromatic_atoms = set([
        atom.GetIdx() for atom in product_mol.GetAtoms()
        if atom.GetIsAromatic()
    ])

    # Check if the number of aromatic atoms changed significantly
    if abs(len(reactant_aromatic_atoms) - len(product_aromatic_atoms)) > 2:
        results["detected_issues"].append(
            f"Significant change in aromaticity: Reactant has {len(reactant_aromatic_atoms)} "
            f"aromatic atoms, Product has {len(product_aromatic_atoms)}")

    # Advanced check for substituent position changes on rings
    check_ring_substituent_positions(reactant_mol, product_mol, results)

    # Check for unnecessary bond formations
    reactant_bonds = Counter(
        [bond.GetBondType() for bond in reactant_mol.GetBonds()])
    product_bonds = Counter(
        [bond.GetBondType() for bond in product_mol.GetBonds()])

    if sum(reactant_bonds.values()) < sum(product_bonds.values()):
        results["detected_issues"].append(
            f"Possible unnecessary bonds formed: Reactant has {sum(reactant_bonds.values())} bonds, "
            f"Product has {sum(product_bonds.values())} bonds")

    return results


def check_ring_substituent_positions(reactant_mol, product_mol, results):
    """
    Detect changes in the position of substituents on aromatic rings.
    
    Parameters:
    -----------
    reactant_mol : RDKit Mol
        RDKit molecule object of the reactant
    product_mol : RDKit Mol
        RDKit molecule object of the product
    results : dict
        Results dictionary to update with findings
    """
    # Get all ring systems in both molecules
    reactant_ring_info = identify_ring_systems(reactant_mol)
    product_ring_info = identify_ring_systems(product_mol)

    # If ring counts mismatch, this is already caught in the main function
    if len(reactant_ring_info) != len(product_ring_info):
        return

    # For each aromatic ring, identify and compare substituent patterns
    for r_idx, reactant_ring in enumerate(reactant_ring_info):
        if not reactant_ring['is_aromatic']:
            continue

        # Find a matching aromatic ring in the product
        matching_rings = [
            p for p in product_ring_info if p['is_aromatic']
            and p['size'] == reactant_ring['size'] and not p['matched']
        ]

        if not matching_rings:
            continue

        product_ring = matching_rings[0]
        product_ring['matched'] = True  # Mark this ring as matched

        # Identify substituents and their positions for both rings
        reactant_substituents = identify_substituents(reactant_mol,
                                                      reactant_ring)
        product_substituents = identify_substituents(product_mol, product_ring)

        # Create signature of each substituent
        reactant_sig = {}
        product_sig = {}

        for subst in reactant_substituents:
            sig = get_substituent_signature(reactant_mol, subst)
            if sig not in reactant_sig:
                reactant_sig[sig] = []
            reactant_sig[sig].append(subst)

        for subst in product_substituents:
            sig = get_substituent_signature(product_mol, subst)
            if sig not in product_sig:
                product_sig[sig] = []
            product_sig[sig].append(subst)

        # Check for position changes of similar substituents
        for sig in set(reactant_sig.keys()).intersection(
                set(product_sig.keys())):
            r_positions = [pos_map[s['position']] for s in reactant_sig[sig]]
            p_positions = [pos_map[s['position']] for s in product_sig[sig]]

            # Sort positions for easier comparison
            r_positions.sort()
            p_positions.sort()

            if r_positions != p_positions:
                # We found a substituent that has changed position
                subst_name = get_friendly_substituent_name(sig)
                results["detected_issues"].append(
                    f"Substituent position change detected: {subst_name} moved from "
                    f"{', '.join(r_positions)} to {', '.join(p_positions)} position(s)"
                )
                results["substituent_position_changes"].append({
                    "substituent":
                    subst_name,
                    "from_positions":
                    r_positions,
                    "to_positions":
                    p_positions
                })


def identify_ring_systems(mol):
    """
    Identify all ring systems in a molecule and their properties.
    
    Parameters:
    -----------
    mol : RDKit Mol
        RDKit molecule object
        
    Returns:
    --------
    list
        List of dictionaries containing ring information
    """
    rings = []
    ring_info = Chem.GetSSSR(mol)

    for idx, ring in enumerate(ring_info):
        ring_atoms = list(ring)
        is_aromatic = all(
            mol.GetAtomWithIdx(atom_idx).GetIsAromatic()
            for atom_idx in ring_atoms)

        rings.append({
            'id': idx,
            'atoms': ring_atoms,
            'size': len(ring_atoms),
            'is_aromatic': is_aromatic,
            'matched':
            False  # Used later for matching rings between reactant and product
        })

    return rings


def identify_substituents(mol, ring_info):
    """
    Identify all substituents attached to a ring and their positions.
    
    Parameters:
    -----------
    mol : RDKit Mol
        RDKit molecule object
    ring_info : dict
        Dictionary containing ring information
        
    Returns:
    --------
    list
        List of dictionaries containing substituent information
    """
    substituents = []
    ring_atoms = set(ring_info['atoms'])

    # Get connections from ring atoms to non-ring atoms
    for ring_atom_idx in ring_atoms:
        ring_atom = mol.GetAtomWithIdx(ring_atom_idx)

        for neighbor in ring_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()

            # Skip atoms that are part of the ring
            if neighbor_idx in ring_atoms:
                continue

            # Determine the position (ortho, meta, para) relative to other substituents
            position = determine_ring_position(mol, ring_atom_idx, ring_atoms,
                                               ring_info['size'])

            # Find the entire substituent group connected to this point
            subst_atoms = get_connected_atoms(mol, neighbor_idx, ring_atoms)

            substituents.append({
                'attachment_point': ring_atom_idx,
                'first_atom': neighbor_idx,
                'atoms': subst_atoms,
                'position': position
            })

    return substituents


def determine_ring_position(mol, atom_idx, ring_atoms, ring_size):
    """
    Determine the position of a substituent on a ring (ortho, meta, para, etc.).
    
    Parameters:
    -----------
    mol : RDKit Mol
        RDKit molecule object
    atom_idx : int
        Index of the ring atom where the substituent is attached
    ring_atoms : set
        Set of atom indices that form the ring
    ring_size : int
        Size of the ring
        
    Returns:
    --------
    str
        Position description ("ortho", "meta", "para", or numbered position)
    """
    # For 6-membered rings, use ortho/meta/para nomenclature
    if ring_size == 6:
        # Find other substituents on the ring
        other_subst = []
        for ring_atom in ring_atoms:
            if ring_atom == atom_idx:
                continue

            atom = mol.GetAtomWithIdx(ring_atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms:
                    other_subst.append(ring_atom)
                    break

        # If no other substituents, just return position number
        if not other_subst:
            return "1"

        # Calculate distance to other substituents
        distances = {}
        for other in other_subst:
            # Use shortest path through the ring
            path = rdmolops.GetShortestPath(mol, atom_idx, other)
            if path:
                path_len = len(
                    path
                ) - 1  # Subtract 1 because path includes both endpoints

                # Convert distance to position name
                if path_len == 1:
                    pos = "ortho"
                elif path_len == 2:
                    pos = "meta"
                elif path_len == 3:
                    pos = "para"
                else:
                    pos = str(path_len)

                distances[other] = pos

        # Return the closest position if multiple are found
        if distances:
            positions = list(distances.values())
            # Prioritize ortho, then meta, then para for consistent naming
            if "ortho" in positions:
                return "ortho"
            elif "meta" in positions:
                return "meta"
            elif "para" in positions:
                return "para"
            else:
                return positions[0]

    # For other ring sizes, use numbered positions (1, 2, 3, etc.)
    return "1"  # Default for now


def get_connected_atoms(mol, start_idx, exclude_atoms):
    """
    Get all atoms connected to a starting atom, excluding a set of atoms.
    
    Parameters:
    -----------
    mol : RDKit Mol
        RDKit molecule object
    start_idx : int
        Index of the starting atom
    exclude_atoms : set
        Set of atom indices to exclude
        
    Returns:
    --------
    list
        List of atom indices that form the connected component
    """
    visited = set([start_idx])
    queue = [start_idx]

    while queue:
        current = queue.pop(0)
        atom = mol.GetAtomWithIdx(current)

        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()

            if neighbor_idx not in visited and neighbor_idx not in exclude_atoms:
                visited.add(neighbor_idx)
                queue.append(neighbor_idx)

    return list(visited)


def get_substituent_signature(mol, substituent):
    """
    Generate a signature for a substituent to identify similar groups.
    
    Parameters:
    -----------
    mol : RDKit Mol
        RDKit molecule object
    substituent : dict
        Dictionary containing substituent information
        
    Returns:
    --------
    str
        Signature string for the substituent
    """
    # Create a fragment of just the substituent
    atoms = substituent['atoms']
    if not atoms:
        return ""

    # Get the SMILES of the fragment
    # This is a simplified approach - a more robust one would create a proper fragment
    atom_symbols = []
    for atom_idx in atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        atom_symbols.append(atom.GetSymbol())

    # Count elements as a basic signature
    elem_counts = Counter(atom_symbols)
    signature = ".".join(f"{elem}{count}"
                         for elem, count in sorted(elem_counts.items()))

    # For more complex substituents, we could use a more sophisticated approach
    # like a Morgan fingerprint or a proper SMILES fragment

    return signature


def get_friendly_substituent_name(signature):
    """
    Convert a substituent signature to a friendly name when possible.
    
    Parameters:
    -----------
    signature : str
        Signature string for the substituent
        
    Returns:
    --------
    str
        Friendly name for the substituent
    """
    # Map of common substituent signatures to friendly names
    common_substituents = {
        "C1": "Methyl",
        "C2": "Ethyl",
        "C3": "Propyl",
        "N1": "Amino",
        "O1": "Hydroxy",
        "O2": "Carboxyl",
        "O2.C1": "Carboxyl acid",
        "Cl1": "Chloro",
        "Br1": "Bromo",
        "F1": "Fluoro",
        "I1": "Iodo",
        "N1.C1": "Methylamino",
        "C1.O1": "Hydroxy methyl",
        "C1.N1": "Aminomethyl",
        "N1.O1": "Nitro",
        "N1.O2": "Nitro",
        "S1": "Thiol"
    }

    return common_substituents.get(signature, f"Group ({signature})")


# Position mapping for consistent naming
pos_map = {
    "1": "position 1",
    "2": "position 2",
    "3": "position 3",
    "4": "position 4",
    "5": "position 5",
    "6": "position 6",
    "ortho": "ortho",
    "meta": "meta",
    "para": "para"
}

# Test with the provided example
if __name__ == "__main__":
    test_reactant = "c1c(CCNC)cc(C(=O)O)cc1"
    test_product = "c1cc(CCN(C)CC)c(C(=O)O)cc1"

    # test_reactant = "C2CC1CCCC1C2"
    # test_product = "C2CCC1CCCCC1C2"

    print("Testing with provided example:")
    print(f"Reactant: {test_reactant}")
    print(f"Product: {test_product}")

    results = compare_molecules(test_reactant, test_product)

    print("\nValidation Results:")
    print("-------------------")

    if results['substituent_position_changes']:
        print("\nSubstituent Position Changes:")
        for change in results['substituent_position_changes']:
            print(
                f"- {change['substituent']} moved from {', '.join(change['from_positions'])} "
                f"to {', '.join(change['to_positions'])}")

    if results['detected_issues']:
        print("\nDetected Issues:")
        for issue in results['detected_issues']:
            print(f"- {issue}")
