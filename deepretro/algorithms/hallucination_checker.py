"""Heuristic hallucination checker for retrosynthetic reaction steps.

When a retrosynthetic step (product → reactant) is proposed, the
predicted reactant may contain structural mistakes: atoms appearing or vanishing, rings changing size,
substituents jumping to a different position on an aromatic ring, and so on.

This module catches those mistakes automatically by comparing the
reactant and product. Two main entry points are provided:

* `hallucination_compare_molecules` — runs every check and returns a
  detailed breakdown of what (if anything) looks wrong.
* `calculate_hallucination_score` — distils the breakdown into a single
  0–100 score (100 = looks fine, 0 = almost certainly hallucinated) with
  a severity label (low / medium / high / critical).
"""

from collections import Counter
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdmolops


# Position mapping for consistent naming
pos_map: dict[str, str] = {
    "1": "position 1",
    "2": "position 2",
    "3": "position 3",
    "4": "position 4",
    "5": "position 5",
    "6": "position 6",
    "ortho": "ortho",
    "meta": "meta",
    "para": "para",
}


def hallucination_compare_molecules(
    reactant_smiles: str,
    product_smiles: str,
) -> dict[str, Any]:
    """Compare a reactant and product molecule to detect potential hallucinations.

    Given two SMILES strings, this function parses both molecules and checks atom-count consistency, ring-size changes, substituent position
    swaps, aromaticity shifts, and unnecessary bond formations.

    Parameters
    ----------
    reactant_smiles : str
        SMILES string of the reactant molecule.
    product_smiles : str
        SMILES string of the product molecule.

    Returns
    -------
    results : dict[str, Any]
        * ``valid_reactant`` (bool) — reactant SMILES parsed OK.
        * ``valid_product`` (bool) — product SMILES parsed OK.
        * ``atom_count_consistent`` (bool) — all elements match.
        * ``ring_size_changes`` (list[str]) — rings added/removed.
        * ``substituent_position_changes`` (list[dict]) — position swaps.
        * ``detected_issues`` (list[str]) — all issues found (empty if clean).

    Examples
    --------
    >>> from deepretro.algorithms import hallucination_compare_molecules
    >>> res = hallucination_compare_molecules("c1ccccc1", "c1ccccc1OC")
    >>> res["valid_reactant"] and res["valid_product"]
    True
    """
    results = {
        "valid_reactant": False,
        "valid_product": False,
        "atom_count_consistent": False,
        "ring_size_changes": [],
        "substituent_position_changes": [],
        "detected_issues": [],
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
    reactant_atoms = Counter([atom.GetSymbol() for atom in reactant_mol.GetAtoms()])
    product_atoms = Counter([atom.GetSymbol() for atom in product_mol.GetAtoms()])

    # Check atom count consistency
    for atom_symbol in set(list(reactant_atoms.keys()) + list(product_atoms.keys())):
        if reactant_atoms.get(atom_symbol, 0) != product_atoms.get(atom_symbol, 0):
            results["detected_issues"].append(
                f"Atom count mismatch for {atom_symbol}: "
                f"Reactant has {reactant_atoms.get(atom_symbol, 0)}, "
                f"Product has {product_atoms.get(atom_symbol, 0)}"
            )

    if not any("Atom count mismatch" in issue for issue in results["detected_issues"]):
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
            f"Product rings {product_ring_sizes}"
        )

        # Report specific ring changes
        for r_size in reactant_ring_sizes:
            if reactant_ring_sizes.count(r_size) > product_ring_sizes.count(r_size):
                results["ring_size_changes"].append(f"{r_size}-membered ring removed")

        for p_size in product_ring_sizes:
            if product_ring_sizes.count(p_size) > reactant_ring_sizes.count(p_size):
                results["ring_size_changes"].append(f"{p_size}-membered ring added")

    # Check for aromatic ring changes
    reactant_aromatic_atoms = set(
        [atom.GetIdx() for atom in reactant_mol.GetAtoms() if atom.GetIsAromatic()]
    )
    product_aromatic_atoms = set(
        [atom.GetIdx() for atom in product_mol.GetAtoms() if atom.GetIsAromatic()]
    )

    # Check if the number of aromatic atoms changed significantly
    if abs(len(reactant_aromatic_atoms) - len(product_aromatic_atoms)) > 2:
        results["detected_issues"].append(
            f"Significant change in aromaticity: Reactant has {len(reactant_aromatic_atoms)} "
            f"aromatic atoms, Product has {len(product_aromatic_atoms)}"
        )

    # Advanced check for substituent position changes on rings
    check_ring_substituent_positions(reactant_mol, product_mol, results)

    # Check for unnecessary bond formations
    reactant_bonds = Counter([bond.GetBondType() for bond in reactant_mol.GetBonds()])
    product_bonds = Counter([bond.GetBondType() for bond in product_mol.GetBonds()])

    if sum(reactant_bonds.values()) < sum(product_bonds.values()):
        results["detected_issues"].append(
            f"Possible unnecessary bonds formed: Reactant has {sum(reactant_bonds.values())} bonds, "
            f"Product has {sum(product_bonds.values())} bonds"
        )

    return results


def check_ring_substituent_positions(
    reactant_mol: Chem.Mol,
    product_mol: Chem.Mol,
    results: dict[str, Any],
) -> None:
    """Detect changes in the position of substituents on aromatic rings.

    For each aromatic ring that appears in both the reactant and the
    product, this function figures out what groups are attached and
    where (ortho / meta / para).  If the same group shows up at a
    different position in the product, that is flagged, it almost
    always means the LLM hallucinated the position.

    Findings are written directly into *results*.

    Parameters
    ----------
    reactant_mol : rdkit.Chem.Mol
        RDKit molecule object of the reactant.
    product_mol : rdkit.Chem.Mol
        RDKit molecule object of the product.
    results : dict
        Results dictionary to update with findings.

    Examples
    --------
    >>> from rdkit import Chem
    >>> r_mol = Chem.MolFromSmiles("c1ccc(O)cc1")   # phenol
    >>> p_mol = Chem.MolFromSmiles("c1ccc(O)cc1")   # same phenol
    >>> res = {"detected_issues": [], "substituent_position_changes": []}
    >>> check_ring_substituent_positions(r_mol, p_mol, res)
    >>> res["substituent_position_changes"]
    []
    """
    # Get all ring systems in both molecules
    reactant_ring_info = identify_ring_systems(reactant_mol)
    product_ring_info = identify_ring_systems(product_mol)

    # If ring counts mismatch, this is already caught in the main function
    if len(reactant_ring_info) != len(product_ring_info):
        return

    # For each aromatic ring, identify and compare substituent patterns
    for r_idx, reactant_ring in enumerate(reactant_ring_info):
        if not reactant_ring["is_aromatic"]:
            continue

        # Find a matching aromatic ring in the product
        matching_rings = [
            p
            for p in product_ring_info
            if p["is_aromatic"]
            and p["size"] == reactant_ring["size"]
            and not p["matched"]
        ]

        if not matching_rings:
            continue

        product_ring = matching_rings[0]
        product_ring["matched"] = True  # Mark this ring as matched

        # Identify substituents and their positions for both rings
        reactant_substituents = identify_substituents(reactant_mol, reactant_ring)
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
        for sig in set(reactant_sig.keys()).intersection(set(product_sig.keys())):
            r_positions = [pos_map[s["position"]] for s in reactant_sig[sig]]
            p_positions = [pos_map[s["position"]] for s in product_sig[sig]]

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
                results["substituent_position_changes"].append(
                    {
                        "substituent": subst_name,
                        "from_positions": r_positions,
                        "to_positions": p_positions,
                    }
                )


def identify_ring_systems(mol: Chem.Mol) -> list[dict[str, Any]]:
    """Identify all ring systems in a molecule and their properties.

    Walks the SSSR (Smallest Set of Smallest Rings) that RDKit computes
    and, for each ring, notes how many atoms it has, which atom indices
    belong to it, and whether every atom in the ring is aromatic.  The
    ``matched`` flag starts as ``False`` and is used later when pairing
    up rings between reactant and product.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.

    Returns
    -------
    rings : list of dict
        Each dict has keys ``id``, ``atoms``, ``size``, ``is_aromatic``,
        and ``matched``.

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("c1ccccc1")
    >>> rings = identify_ring_systems(mol)
    >>> len(rings)
    1
    >>> rings[0]["size"]
    6
    >>> rings[0]["is_aromatic"]
    True
    """
    rings = []
    ring_info = Chem.GetSSSR(mol)

    for idx, ring in enumerate(ring_info):
        ring_atoms = list(ring)
        is_aromatic = all(
            mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring_atoms
        )

        rings.append(
            {
                "id": idx,
                "atoms": ring_atoms,
                "size": len(ring_atoms),
                "is_aromatic": is_aromatic,
                "matched": False,  # Used later for matching rings between reactant and product
            }
        )

    return rings


def identify_substituents(
    mol: Chem.Mol,
    ring_info: dict[str, Any],
) -> list[dict[str, Any]]:
    """Identify all substituents attached to a ring and their positions.

    Walks the atoms of the ring and, for every neighbour that is *not*
    part of the ring, traces out the full substituent group and labels
    its attachment point as ortho / meta / para (for 6-membered rings)
    or a numbered position (for other ring sizes).

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    ring_info : dict
        Ring descriptor as returned by `identify_ring_systems`.

    Returns
    -------
    substituents : list of dict
        Each dict has keys ``attachment_point``, ``first_atom``,
        ``atoms``, and ``position``.

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("c1ccc(O)cc1")  # phenol
    >>> rings = identify_ring_systems(mol)
    >>> subs = identify_substituents(mol, rings[0])
    >>> len(subs) >= 1
    True
    >>> subs[0]["position"] in ("1", "ortho", "meta", "para")
    True
    """
    substituents = []
    ring_atoms = set(ring_info["atoms"])

    # Get connections from ring atoms to non-ring atoms
    for ring_atom_idx in ring_atoms:
        ring_atom = mol.GetAtomWithIdx(ring_atom_idx)

        for neighbor in ring_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()

            # Skip atoms that are part of the ring
            if neighbor_idx in ring_atoms:
                continue

            # Determine the position (ortho, meta, para) relative to other substituents
            position = determine_ring_position(
                mol, ring_atom_idx, ring_atoms, ring_info["size"]
            )

            # Find the entire substituent group connected to this point
            subst_atoms = get_connected_atoms(mol, neighbor_idx, ring_atoms)

            substituents.append(
                {
                    "attachment_point": ring_atom_idx,
                    "first_atom": neighbor_idx,
                    "atoms": subst_atoms,
                    "position": position,
                }
            )

    return substituents


def determine_ring_position(
    mol: Chem.Mol,
    atom_idx: int,
    ring_atoms: set[int],
    ring_size: int,
) -> str:
    """
    Determine the position of a substituent on a ring.

    For 6-membered rings uses ortho/meta/para nomenclature.
    For other ring sizes returns numbered positions.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    atom_idx : int
        Index of the ring atom the substituent is bonded to.
    ring_atoms : set[int]
        All atom indices that belong to the ring.
    ring_size : int
        Size of the ring.

    Returns
    -------
    position : str
        ``"ortho"``, ``"meta"``, ``"para"``, or a numbered position.

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("c1ccc(O)cc1")  # phenol
    >>> ring_atoms = set(range(6))
    >>> pos = determine_ring_position(mol, 3, ring_atoms, 6)
    >>> pos in ("1", "ortho", "meta", "para")
    True
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
                path_len = (
                    len(path) - 1
                )  # Subtract 1 because path includes both endpoints

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


def get_connected_atoms(
    mol: Chem.Mol,
    start_idx: int,
    exclude_atoms: set[int],
) -> list[int]:
    """
    Get all atoms connected to a starting atom, excluding a set of atoms.

    Starting from *start_idx* (typically the first atom outside a ring),
    this does a breadth-first walk along bonds and collects every atom
    it reaches. It will *not* cross into any atom listed in
    *exclude_atoms*, this is how we stop at the ring boundary and only
    get the substituent itself.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    start_idx : int
        Atom index to start the walk from.
    exclude_atoms : set[int]
        Atom indices to treat as barriers (usually the ring atoms).

    Returns
    -------
    atoms : list of int
        List of atom indices that form the connected component.

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("c1ccc(OC)cc1")  # methoxybenzene
    >>> ring_atoms = set(range(6))
    >>> # atom 6 is the O attached to the ring; BFS from there excluding ring
    >>> connected = get_connected_atoms(mol, 6, ring_atoms)
    >>> len(connected) >= 1
    True
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


def get_substituent_signature(
    mol: Chem.Mol,
    substituent: dict[str, Any],
) -> str:
    """
    Generate a signature for a substituent to identify similar groups.

    Counts element types in the substituent atoms and returns a sorted
    dot-separated string (e.g. ``"C2.O1"``).

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    substituent : dict
        Substituent descriptor (must contain an ``atoms`` key with
        a list of atom indices).

    Returns
    -------
    signature : str
        Signature string for the substituent.

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("c1ccc(O)cc1")  # phenol
    >>> subst = {"atoms": [6]}  # the oxygen atom
    >>> get_substituent_signature(mol, subst)
    'O1'
    """
    # Create a fragment of just the substituent
    atoms = substituent["atoms"]
    if not atoms:
        return ""

    # Get the SMILES of the fragment
    # This is a simplified approach, a more robust one would create a proper fragment
    atom_symbols = []
    for atom_idx in atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        atom_symbols.append(atom.GetSymbol())

    # Count elements as a basic signature
    elem_counts = Counter(atom_symbols)
    signature = ".".join(
        f"{elem}{count}" for elem, count in sorted(elem_counts.items())
    )

    # For more complex substituents, we could use a more sophisticated approach
    # like a Morgan fingerprint or a proper SMILES fragment

    return signature


def get_friendly_substituent_name(signature: str) -> str:
    """
    Convert a substituent signature to a friendly name when possible.

    Parameters
    ----------
    signature : str
        Element-count signature (e.g. ``"C1"``, ``"N1.O2"``).

    Returns
    -------
    name : str
        Friendly name (e.g. ``"Methyl"``), or ``"Group (<signature>)"``
        if no match is found.

    Examples
    --------
    >>> get_friendly_substituent_name("C1")
    'Methyl'
    >>> get_friendly_substituent_name("Br1")
    'Bromo'
    >>> get_friendly_substituent_name("X99")
    'Group (X99)'
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
        "S1": "Thiol",
    }

    return common_substituents.get(signature, f"Group ({signature})")


def calculate_hallucination_score(
    reactant_smiles: str,
    product_smiles: str,
) -> dict[str, Any]:
    """
    Calculate a hallucination score for a chemical transformation.

    This is the high-level entry point.  It runs
    `hallucination_compare_molecules` under the hood and then converts
    each kind of issue into a point deduction from a perfect score of 100.
    Bigger problems cost more points (e.g. a substituent jumping position
    costs 60, while one extra bond costs only 5).  The final score is
    clamped to 0–100 and labelled with a severity:

    * **≥ 80** → ``"low"`` — looks plausible
    * **40–79** → ``"medium"`` — worth a second look
    * **20–39** → ``"high"`` — likely hallucinated
    * **< 20** → ``"critical"`` — almost certainly wrong

    If either SMILES string cannot be parsed, the score is 0 / critical.

    Parameters
    ----------
    reactant_smiles : str
        SMILES string of the reactant molecule.
    product_smiles : str
        SMILES string of the product molecule.

    Returns
    -------
    result : dict
        Dictionary with keys ``score`` (int, 0–100), ``severity``
        (``"low"`` / ``"medium"`` / ``"high"`` / ``"critical"``),
        ``penalties`` (list of str), and ``message`` (str).

    Examples
    --------
    >>> from deepretro.algorithms import calculate_hallucination_score
    >>> result = calculate_hallucination_score("c1ccccc1", "c1ccccc1OC")
    >>> result["severity"]
    'low'
    >>> result["score"] >= 80
    True
    """
    # Get the detailed comparison results first
    comparison_results = hallucination_compare_molecules(
        reactant_smiles, product_smiles
    )

    # Initialize the score at 100 (no hallucinations)
    base_score = 100
    penalty_factors = []
    penalty_descriptions = []

    # Check if molecules are valid - severe penalty if not
    if (
        not comparison_results["valid_reactant"]
        or not comparison_results["valid_product"]
    ):
        return {
            "score": 0,
            "severity": "critical",
            "message": "Invalid SMILES string detected - cannot assess transformation",
        }

    # Apply penalties based on detected issues

    # 1. Atom count consistency - Critical issue
    if not comparison_results["atom_count_consistent"]:
        atom_mismatch_penalties = []

        for issue in comparison_results["detected_issues"]:
            if "Atom count mismatch" in issue:
                # Extract the difference in atom counts
                parts = issue.split(":")[1].strip()
                reactant_count = int(parts.split(",")[0].split()[-1])
                product_count = int(parts.split(",")[1].split()[-1])
                difference = abs(reactant_count - product_count)

                # Penalty: 5 points per atom mismatch
                penalty = min(5 * difference, 100)
                atom_mismatch_penalties.append(penalty)

                penalty_descriptions.append(
                    f"Atom count inconsistency: -{penalty} points"
                )

        # Take the maximum penalty from atom mismatches
        if atom_mismatch_penalties:
            penalty_factors.append(max(atom_mismatch_penalties))

    # 2. Ring size changes - Potential issue, but could be valid in some reactions
    if comparison_results["ring_size_changes"]:
        # Check how many ring changes occurred
        num_ring_changes = len(comparison_results["ring_size_changes"])

        # Penalty: 25 points per ring change
        ring_penalty = min(25 * num_ring_changes, 50)
        penalty_factors.append(ring_penalty)
        penalty_descriptions.append(f"Ring structure changes: -{ring_penalty} points")

    # 3. Substituent position changes - Usually suspicious
    if comparison_results["substituent_position_changes"]:
        # Check how many substituent position changes
        num_position_changes = len(comparison_results["substituent_position_changes"])

        # Penalty: 60 points per substituent position change
        position_penalty = min(60 * num_position_changes, 100)
        penalty_factors.append(position_penalty)
        penalty_descriptions.append(
            f"Substituent position changes: -{position_penalty} points"
        )

    # 4. Aromaticity changes - Significant structural change
    for issue in comparison_results["detected_issues"]:
        if "Significant change in aromaticity" in issue:
            aromaticity_penalty = 40
            penalty_factors.append(aromaticity_penalty)
            penalty_descriptions.append(
                f"Significant aromaticity changes: -{aromaticity_penalty} points"
            )

    # 5. Unnecessary bond formations - Could indicate hallucination
    for issue in comparison_results["detected_issues"]:
        if "Possible unnecessary bonds formed" in issue:
            # Extract number of additional bonds
            parts = issue.split(":")[1].strip()
            reactant_bonds = int(parts.split(",")[0].split()[-2])
            product_bonds = int(parts.split(",")[1].split()[-2])
            additional_bonds = product_bonds - reactant_bonds

            # Penalty: 5 points per additional bond
            bond_penalty = min(5 * additional_bonds, 30)
            penalty_factors.append(bond_penalty)
            penalty_descriptions.append(
                f"Unnecessary bond formations: -{bond_penalty} points"
            )

    # Calculate the final score by applying all penalties
    final_score = base_score
    for penalty in penalty_factors:
        final_score -= penalty

    # Ensure score doesn't go below 0
    final_score = max(0, final_score)

    # Determine severity level based on score
    if final_score >= 80:
        severity = "low"
    elif final_score >= 55:
        severity = "low-medium"
    elif final_score >= 40:
        severity = "medium"
    elif final_score >= 20:
        severity = "high"
    else:
        severity = "critical"

    return {
        "score": final_score,
        "severity": severity,
        "penalties": penalty_descriptions,
        "message": interpret_score(final_score),
    }


def interpret_score(score: int) -> str:
    """Turn a numeric hallucination score into a sentence a non-expert can read.

    This is called automatically by `calculate_hallucination_score` to
    fill the ``message`` field, but you can also use it standalone if
    you already have a score.

    Parameters
    ----------
    score : int
        Hallucination score (0 = worst, 100 = best).

    Returns
    -------
    message : str
        One-sentence plain-English interpretation.

    Examples
    --------
    >>> from deepretro.algorithms import interpret_score
    >>> interpret_score(95)
    'Highly reliable transformation with minimal or no structural inconsistencies'
    """
    if score >= 90:
        return "Highly reliable transformation with minimal or no structural inconsistencies"
    elif score >= 80:
        return "Generally reliable transformation with minor structural inconsistencies"
    elif score >= 70:
        return "Mostly reliable transformation with some structural inconsistencies"
    elif score >= 50:
        return "Questionable transformation with significant structural inconsistencies"
    elif score >= 30:
        return "Likely hallucination with major structural inconsistencies"
    elif score >= 10:
        return "Severe hallucination with critical structural inconsistencies"
    else:
        return "Complete hallucination or invalid transformation"
