"""Structural consistency heuristics for retrosynthetic reaction steps.

Compare proposed reactants with a product, then score the resulting atom, ring,
substituent and bond differences. Scores range from 0 to 100 and are not
calibrated probabilities of chemical correctness. Reaction conditions, mechanism,
stereoselectivity and experimental feasibility require separate assessment.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdmolops

from deepretro.algorithms.hallucination_weights import (
    DEFAULT_WEIGHTS,
    HallucinationWeights,
)

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


def _position_label(raw: str) -> str:
    """Human-readable form of a raw position token.

    ``pos_map`` only knows the six numbered slots and ortho/meta/para. Ring-bond
    separations from non-six-membered rings are formatted here instead of
    silently missing from the map.

    Examples
    --------
    >>> _position_label("ortho")
    'ortho'
    >>> _position_label("sep2")
    'separation 2'
    """
    if raw.startswith("sep"):
        return f"separation {raw[3:]}"
    return pos_map.get(raw, raw)


#: Minimum heavy-atom overlap for retaining a potential contributing fragment.
#: This is a heuristic threshold, not a validated assignment of reagent roles.
REAGENT_MCS_MIN_FRACTION = 0.5


def _normalise_isotopes(mol: Chem.Mol) -> Chem.Mol:
    """Strip isotope labels and explicit hydrogens from a copy of *mol*.

    An isotopic label does not change which skeleton was proposed, but RDKit
    represents ``[2H]`` as an explicit atom and ``[13C]`` as an ordinary carbon,
    so the two were scored inconsistently.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to normalise. Not modified in place.

    Returns
    -------
    rdkit.Chem.Mol
        The normalised molecule, or *mol* unchanged if normalisation fails.

    Examples
    --------
    >>> callable(_normalise_isotopes)
    True
    """
    try:
        copy = Chem.Mol(mol)
        for atom in copy.GetAtoms():
            atom.SetIsotope(0)
        return Chem.RemoveHs(copy)
    except Exception:  # noqa: BLE001 - normalisation must never fail the check
        return mol


def _drop_spectator_reagents(reactant_smiles: str, product_mol: Chem.Mol | None) -> str:
    """Heuristically omit fragments with little structural overlap to the product.

    Exact-bond maximum common substructures identify possible spectators. Keep
    uncertain fragments, and retain the complete input when the remaining heavy
    atoms cannot cover the product. This does not establish atom provenance.

    Parameters
    ----------
    reactant_smiles : str
        Dot-joined reactant SMILES.
    product_mol : rdkit.Chem.Mol
        Parsed product.

    Returns
    -------
    str
        The contributing fragments, dot-joined. The input is returned unchanged
        when it has one fragment, when the product will not parse, or when every
        fragment looks like a spectator.

    Examples
    --------
    >>> callable(_drop_spectator_reagents)
    True
    """
    fragments = reactant_smiles.split(".")
    if len(fragments) < 2 or product_mol is None:
        return reactant_smiles

    keep: list[str] = []
    for fragment in fragments:
        mol = Chem.MolFromSmiles(fragment)
        if mol is None:
            keep.append(fragment)  # unparseable is the validity check's business
            continue
        if mol.GetNumHeavyAtoms() == 0:
            continue
        try:
            result = rdFMCS.FindMCS(
                [mol, product_mol],
                timeout=1,
                ringMatchesRingOnly=True,
                # Without exact bond matching an aliphatic ring matches an
                # aromatic one, and DCC's cyclohexyl groups "contribute" to a
                # benzyl product.
                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
            )
            if result.canceled:
                # A timeout provides no evidence that a fragment is a spectator.
                covered = 1.0
            elif result.numAtoms == 0:
                covered = 0.0
            else:
                covered = result.numAtoms / mol.GetNumHeavyAtoms()
        except Exception:  # noqa: BLE001 - never fail the check over MCS trouble
            covered = 1.0
        if covered >= REAGENT_MCS_MIN_FRACTION:
            keep.append(fragment)

    # Everything looked like a spectator: keep the input rather than compare the
    # product against nothing.
    if not keep:
        return reactant_smiles

    kept_atoms = 0
    for fragment in keep:
        mol = Chem.MolFromSmiles(fragment)
        if mol is not None:
            kept_atoms += mol.GetNumHeavyAtoms()
    if kept_atoms < product_mol.GetNumHeavyAtoms():
        return reactant_smiles

    return ".".join(keep)


def hallucination_compare_molecules(
    reactant_smiles: str,
    product_smiles: str,
    *,
    arom_trigger: int = 2,
    strip_reagents: bool = True,
) -> dict[str, Any]:
    """Compare a reactant and product molecule to detect potential hallucinations.

    Parse both molecules and compare atom counts, ring sizes, substituent
    positions, aromaticity and bond counts. Findings are structural warnings.

    Parameters
    ----------
    reactant_smiles : str
        SMILES string of the reactant molecule.
    product_smiles : str
        SMILES string of the product molecule.
    arom_trigger : int, optional
        Aromatic-atom count difference above which to report an issue.
    strip_reagents : bool, optional
        Apply the conservative spectator heuristic before skeleton comparison.

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
        # Numeric fields -- the only thing score_from_comparison reads.
        "atom_count_deltas": {},
        "aromatic_atom_delta": 0,
        "bond_delta": 0,
        "n_ring_changes": 0,
        "n_position_changes": 0,
        # Initialised here so the dict has one shape on every path. These are
        # written to feature files and JSONL, where a column that appears and
        # disappears depending on the input is a real nuisance.
        "degenerate": False,
        "carbon_radical": False,
        "dropped_reagents": False,
    }

    # Check if SMILES strings are valid
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    if reactant_mol is None or reactant_mol.GetNumAtoms() == 0:
        results["detected_issues"].append("Invalid reactant SMILES string")
        return results
    else:
        results["valid_reactant"] = True

    if product_mol is None or product_mol.GetNumAtoms() == 0:
        results["detected_issues"].append("Invalid product SMILES string")
        return results
    else:
        results["valid_product"] = True

    # Strip spectators before isotope normalization; detect no-ops while
    # retaining isotope identity. Atom maps are correspondence metadata.
    for mol in (reactant_mol, product_mol):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    reactant_mol = Chem.MolFromSmiles(Chem.MolToSmiles(reactant_mol))
    product_mol = Chem.MolFromSmiles(Chem.MolToSmiles(product_mol))
    reactant_smiles = Chem.MolToSmiles(reactant_mol)

    if strip_reagents:
        stripped = _drop_spectator_reagents(reactant_smiles, product_mol)
        if stripped != reactant_smiles:
            reduced = Chem.MolFromSmiles(stripped)
            if reduced is not None:
                reactant_mol = reduced
                results["dropped_reagents"] = True

    # A step whose reactants ARE the product is not a retrosynthesis step: it
    # decomposes nothing. Compared BEFORE isotope normalisation so that a
    # labelling step, whose skeletons match but whose isotopes differ, is not
    # swept up.
    results["degenerate"] = product_mol.GetNumHeavyAtoms() == 0 or (
        Chem.MolToSmiles(reactant_mol) == Chem.MolToSmiles(product_mol)
    )

    # Isotope labels and explicit hydrogens are graph noise, not skeleton
    # changes. RDKit models [2H] as its own atom but [13C] as an ordinary
    # carbon, so d10-ethylbenzene -> ethylbenzene scored 50 while a 13C label
    # scored 100. Normalising both makes the treatment consistent.
    reactant_mol = _normalise_isotopes(reactant_mol)
    product_mol = _normalise_isotopes(product_mol)

    results["carbon_radical"] = any(
        atom.GetSymbol() == "C" and atom.GetNumRadicalElectrons() > 0
        for mol in (reactant_mol, product_mol)
        for atom in mol.GetAtoms()
    )

    # Get basic molecule properties
    reactant_atoms = Counter([atom.GetSymbol() for atom in reactant_mol.GetAtoms()])
    product_atoms = Counter([atom.GetSymbol() for atom in product_mol.GetAtoms()])

    # Check atom count consistency
    for atom_symbol in sorted(set(reactant_atoms) | set(product_atoms)):
        r_count = reactant_atoms.get(atom_symbol, 0)
        p_count = product_atoms.get(atom_symbol, 0)
        if r_count != p_count:
            results["atom_count_deltas"][atom_symbol] = abs(r_count - p_count)
            results["detected_issues"].append(
                f"Atom count mismatch for {atom_symbol}: "
                f"Reactant has {r_count}, "
                f"Product has {p_count}"
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

        reactant_counts = Counter(reactant_ring_sizes)
        product_counts = Counter(product_ring_sizes)
        for size, count in sorted((reactant_counts - product_counts).items()):
            results["ring_size_changes"].extend(
                [f"{size}-membered ring removed"] * count
            )
        for size, count in sorted((product_counts - reactant_counts).items()):
            results["ring_size_changes"].extend([f"{size}-membered ring added"] * count)

    # Check for aromatic ring changes
    reactant_aromatic_atoms = set(
        [atom.GetIdx() for atom in reactant_mol.GetAtoms() if atom.GetIsAromatic()]
    )
    product_aromatic_atoms = set(
        [atom.GetIdx() for atom in product_mol.GetAtoms() if atom.GetIsAromatic()]
    )

    # Check if the number of aromatic atoms changed significantly
    aromatic_delta = abs(len(reactant_aromatic_atoms) - len(product_aromatic_atoms))
    results["aromatic_atom_delta"] = aromatic_delta
    if aromatic_delta > arom_trigger:
        results["detected_issues"].append(
            f"Significant change in aromaticity: Reactant has {len(reactant_aromatic_atoms)} "
            f"aromatic atoms, Product has {len(product_aromatic_atoms)}"
        )

    # Advanced check for substituent position changes on rings
    check_ring_substituent_positions(reactant_mol, product_mol, results)

    # Check for unnecessary bond formations
    reactant_bonds = Counter([bond.GetBondType() for bond in reactant_mol.GetBonds()])
    product_bonds = Counter([bond.GetBondType() for bond in product_mol.GetBonds()])

    r_bonds = sum(reactant_bonds.values())
    p_bonds = sum(product_bonds.values())
    results["bond_delta"] = max(0, p_bonds - r_bonds)
    if r_bonds < p_bonds:
        results["detected_issues"].append(
            f"Possible unnecessary bonds formed: Reactant has {r_bonds} bonds, "
            f"Product has {p_bonds} bonds"
        )

    results["n_ring_changes"] = len(results["ring_size_changes"])
    results["n_position_changes"] = len(results["substituent_position_changes"])

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
    different position in the product, a structural warning is recorded.
    Ring matching is heuristic and does not establish an atom correspondence.

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
    for reactant_ring in reactant_ring_info:
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

        # Pair with the ring whose substituents overlap most, not simply the
        # first same-sized one. With matching_rings[0] the verdict depended on
        # fragment order: the same two reactants scored 75 one way round and 0
        # the other.
        reactant_substituents = identify_substituents(reactant_mol, reactant_ring)
        reactant_sigs = {
            get_substituent_signature(reactant_mol, s) for s in reactant_substituents
        }

        def _overlap(candidate: dict) -> tuple[int, int]:
            subs = identify_substituents(product_mol, candidate)
            sigs = {get_substituent_signature(product_mol, s) for s in subs}
            # Prefer most shared substituents, then the closest count.
            return (
                len(reactant_sigs & sigs),
                -abs(len(subs) - len(reactant_substituents)),
            )

        product_ring = max(matching_rings, key=_overlap)
        product_ring["matched"] = True  # Mark this ring as matched

        product_substituents = identify_substituents(product_mol, product_ring)

        if len(reactant_substituents) != len(product_substituents):
            continue

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
        for sig in sorted(set(reactant_sig) & set(product_sig)):
            r_positions = [_position_label(s["position"]) for s in reactant_sig[sig]]
            p_positions = [_position_label(s["position"]) for s in product_sig[sig]]

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

    other_subst = []
    for ring_atom in ring_atoms:
        if ring_atom == atom_idx:
            continue
        for neighbor in mol.GetAtomWithIdx(ring_atom).GetNeighbors():
            if neighbor.GetIdx() not in ring_atoms:
                other_subst.append(ring_atom)
                break
    if not other_subst:
        return "1"

    ring_set = set(ring_atoms)
    distances = []
    for other in other_subst:
        path = rdmolops.GetShortestPath(mol, atom_idx, other)
        # Only count paths that stay on the ring; a path leaving it is not a
        # ring separation.
        if path and set(path) <= ring_set:
            distances.append(len(path) - 1)
    if not distances:
        return "1"
    return f"sep{min(distances)}"


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
    >>> subst = {"atoms": [4]}  # the oxygen atom in this SMILES
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


def score_from_comparison(
    comparison: dict[str, Any],
    weights: HallucinationWeights = DEFAULT_WEIGHTS,
) -> dict[str, Any]:
    """Turn a comparison dict into a 0-100 score under the given weights.

    Pure arithmetic. Reads only the numeric fields -- never ``detected_issues``
    -- so a caller may compute the comparison once and re-score it under many
    weight vectors without touching RDKit again.

    Parameters
    ----------
    comparison : dict
        Output of :func:`hallucination_compare_molecules`.
    weights : HallucinationWeights, optional
        Defaults to the version-2 configuration.

    Returns
    -------
    dict
        ``score`` (int, 0-100), ``severity``, ``penalties``, ``message``. On an
        unassessable comparison, ``unassessable`` is true and penalties are empty.
        ``penalty_total`` retains the sum before the score is clamped at zero.

    Examples
    --------
    >>> from deepretro.algorithms.hallucination_checker import (
    ...     hallucination_compare_molecules, score_from_comparison)
    >>> cmp_ = hallucination_compare_molecules("c1ccccc1", "c1ccccc1")
    >>> score_from_comparison(cmp_)["score"]
    0
    """
    if not comparison.get("valid_reactant") or not comparison.get("valid_product"):
        return {
            "score": 0,
            "penalty_total": 0,
            "unassessable": True,
            "severity": "critical",
            "penalties": [],
            "message": "Invalid SMILES string detected - cannot assess transformation",
        }

    if comparison.get("degenerate"):
        return {
            "score": 0,
            "penalty_total": 0,
            "unassessable": True,
            "severity": "critical",
            "penalties": [],
            "message": (
                "Degenerate step - the reactants are the product, so nothing is "
                "decomposed"
            ),
        }

    # Treated as a validity failure rather than a penalty: the structure parsed,
    # but a carbon radical means the SMILES was truncated, so the atom and ring
    # counts describe a molecule nobody proposed.
    if comparison.get("carbon_radical"):
        return {
            "score": 0,
            "penalty_total": 0,
            "unassessable": True,
            "severity": "critical",
            "penalties": [],
            "message": (
                "Carbon radical detected - radical chemistry is outside the "
                "scope of this heuristic"
            ),
        }

    penalty_factors: list[int] = []
    penalty_descriptions: list[str] = []

    deltas = comparison.get("atom_count_deltas") or {}
    if deltas:
        per_element = [
            min(weights.w_atom * delta, weights.cap_atom) for delta in deltas.values()
        ]
        applied = min(sum(per_element), weights.cap_atom)
        penalty_factors.append(applied)
        # Report what was APPLIED, not what each element would have cost. The
        # per-element lines summed to more than the score whenever cap_atom
        # bound: a case reporting -25, -30 and -30 had only 30 applied, so the
        # explanation overstated the penalty by 55 points.
        if sum(per_element) > applied:
            penalty_descriptions.append(
                f"Atom count inconsistency across {len(per_element)} element(s), "
                f"capped: -{applied} points"
            )
        else:
            for penalty in per_element:
                penalty_descriptions.append(
                    f"Atom count inconsistency: -{penalty} points"
                )

    # 2. Ring size changes.
    n_ring = comparison.get("n_ring_changes", 0)
    if n_ring:
        ring_penalty = min(weights.w_ring * n_ring, weights.cap_ring)
        penalty_factors.append(ring_penalty)
        penalty_descriptions.append(f"Ring structure changes: -{ring_penalty} points")

    # 3. Substituent position changes.
    n_pos = comparison.get("n_position_changes", 0)
    if n_pos:
        position_penalty = min(weights.w_pos * n_pos, weights.cap_pos)
        penalty_factors.append(position_penalty)
        penalty_descriptions.append(
            f"Substituent position changes: -{position_penalty} points"
        )

    # 4. Aromaticity changes.
    if comparison.get("aromatic_atom_delta", 0) > weights.arom_trigger:
        penalty_factors.append(weights.w_arom)
        penalty_descriptions.append(
            f"Significant aromaticity changes: -{weights.w_arom} points"
        )

    # 5. Unnecessary bond formations.
    bond_delta = comparison.get("bond_delta", 0)
    if bond_delta > 0:
        bond_penalty = min(weights.w_bond * bond_delta, weights.cap_bond)
        penalty_factors.append(bond_penalty)
        penalty_descriptions.append(
            f"Unnecessary bond formations: -{bond_penalty} points"
        )

    final_score = max(0, 100 - sum(penalty_factors))

    if final_score >= weights.cut_low:
        severity = "low"
    elif final_score >= weights.cut_medium:
        severity = "medium"
    elif final_score >= weights.cut_high:
        severity = "high"
    else:
        severity = "critical"

    return {
        "score": final_score,
        "penalty_total": int(sum(penalty_factors)),
        "severity": severity,
        "penalties": penalty_descriptions,
        # Derived from the severity that was actually assigned, so the wording
        # cannot contradict it. interpret_score hardcoded 90/80/70/50/30/10
        # while severity uses the tunable cutoffs, which produced reports like
        # "score 30, severity low, Likely hallucination".
        "message": describe_severity(severity, final_score),
    }


def calculate_hallucination_score(
    reactant_smiles: str,
    product_smiles: str,
    weights: HallucinationWeights | None = None,
) -> dict[str, Any]:
    """Score structural consistency from 0 (most penalized) to 100 (unpenalized).

    Parameters
    ----------
    reactant_smiles : str
        Dot-joined SMILES of the proposed reactants.
    product_smiles : str
        SMILES of the product being disconnected.
    weights : HallucinationWeights, optional
        When ``None`` (the default), use the version-2 configuration.

    Returns
    -------
    dict
        Keys ``score`` (int, 0-100), ``severity``
        (``"low"``/``"medium"``/``"high"``/``"critical"``), ``penalties``
        (list of str), ``penalty_total`` and ``message`` (str). An additional
        ``unassessable=True`` marks invalid, empty, no-op or radical inputs.

    Examples
    --------
    >>> from deepretro.algorithms import calculate_hallucination_score
    >>> result = calculate_hallucination_score("c1ccccc1", "c1ccccc1")
    >>> result["severity"]
    'critical'
    >>> result["unassessable"]
    True
    """
    w = DEFAULT_WEIGHTS if weights is None else weights
    comparison = hallucination_compare_molecules(
        reactant_smiles, product_smiles, arom_trigger=w.arom_trigger
    )
    return score_from_comparison(comparison, w)


def describe_severity(severity: str, score: int) -> str:
    """Describe a verdict using the severity that was actually assigned.

    Parameters
    ----------
    severity : str
        One of ``"low"``, ``"medium"``, ``"high"``, ``"critical"``.
    score : int
        The 0-100 score, quoted in the message.

    Returns
    -------
    str
        Human-readable wording consistent with *severity*.

    Examples
    --------
    >>> describe_severity("low", 100)
    'Transformation looks consistent (score 100/100)'
    """
    wording = {
        "low": "Transformation looks consistent",
        "medium": "Minor inconsistencies; review recommended",
        "high": "Significant structural inconsistencies; review required",
        "critical": "Severe structural inconsistencies; review required",
    }
    return f"{wording.get(severity, 'Unrecognised severity')} (score {score}/100)"


def interpret_score(score: int) -> str:
    """Turn a numeric hallucination score into a sentence a non-expert can read.

    This legacy helper uses fixed display thresholds. Scoring uses
    :func:`describe_severity` with the configured severity cutoffs instead.
    Neither description establishes chemical feasibility.

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
    'Minimal or no structural inconsistencies detected'
    """
    if score >= 90:
        return "Minimal or no structural inconsistencies detected"
    elif score >= 80:
        return "Minor structural inconsistencies detected"
    elif score >= 70:
        return "Some structural inconsistencies detected; review recommended"
    elif score >= 50:
        return "Significant structural inconsistencies detected; review recommended"
    elif score >= 30:
        return "Major structural inconsistencies detected; review required"
    elif score >= 10:
        return "Critical structural inconsistencies detected; review required"
    else:
        return "Fundamental structural inconsistencies or unsupported transformation"
