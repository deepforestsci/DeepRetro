"""Molecule helpers for validation, filtering, and simple descriptors."""

from __future__ import annotations

from collections.abc import Sequence

import structlog
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

logger = structlog.get_logger()

_INVALID_SMILES_ERROR = "Invalid SMILES string provided."


def _parse_molecule(smiles: str, *, message: str | None = None):
    """Parse a SMILES string into an RDKit molecule."""

    try:
        molecule = Chem.MolFromSmiles(smiles)
    except Exception:
        molecule = None

    if molecule is None and message is not None:
        logger.warning(message, smiles=smiles)

    return molecule


def is_valid_smiles(smiles: str) -> bool:
    """Check whether a SMILES string can be parsed successfully.

    Parameters
    ----------
    smiles : str
        SMILES string to validate.

    Returns
    -------
    bool
        ``True`` when the SMILES string parses to an RDKit molecule,
        otherwise ``False``.

    Examples
    --------
    >>> is_valid_smiles("CCO")
    True
    >>> is_valid_smiles("not_a_smiles")
    False
    """

    return _parse_molecule(smiles) is not None


def substructure_matching(target_smiles: str, query_smiles: str) -> int:
    """Check whether a query molecule is a substructure of a target molecule.

    Parameters
    ----------
    target_smiles : str
        SMILES string of the target molecule.
    query_smiles : str
        SMILES string of the query molecule.

    Returns
    -------
    int
        ``1`` if the query is a substructure of the target, otherwise ``0``.

    Examples
    --------
    >>> substructure_matching("CCc1ccccc1", "c1ccccc1")
    1
    >>> substructure_matching("CCO", "c1ccccc1")
    0
    """

    target_molecule = _parse_molecule(
        target_smiles,
        message="Unable to parse target molecule for substructure matching",
    )
    query_molecule = _parse_molecule(
        query_smiles,
        message="Unable to parse query molecule for substructure matching",
    )
    if target_molecule is None or query_molecule is None:
        return 0

    return int(target_molecule.HasSubstructMatch(query_molecule))


def validity_check(
    molecule: str,
    res_molecules: Sequence[Sequence[str] | str],
    res_explanations: Sequence[str],
    res_confidence: Sequence[float],
) -> tuple[list[list[str]], list[str], list[float]]:
    """Filter proposed retrosynthesis pathways down to valid precursor sets.

    Parameters
    ----------
    molecule : str
        Target molecule for retrosynthesis.
    res_molecules : Sequence[Sequence[str] | str]
        Candidate precursor pathways returned by the model.
    res_explanations : Sequence[str]
        Explanation for each candidate pathway.
    res_confidence : Sequence[float]
        Confidence score for each candidate pathway.

    Returns
    -------
    tuple[list[list[str]], list[str], list[float]]
        Valid precursor pathways, explanations, and confidence scores. A
        pathway is kept only when every precursor is valid, is not identical
        to the target molecule, and is not a substructure of the target.

    Examples
    --------
    >>> original_logger = validity_check.__globals__["logger"]
    >>> class _SilentLogger:
    ...     def info(self, *args, **kwargs):
    ...         pass
    ...
    ...     def warning(self, *args, **kwargs):
    ...         pass
    >>> validity_check.__globals__["logger"] = _SilentLogger()
    >>> validity_check(
    ...     molecule="c1ccccc1",
    ...     res_molecules=[["CCO", "CCCl"]],
    ...     res_explanations=["valid pathway"],
    ...     res_confidence=[0.9],
    ... )
    ([['CCO', 'CCCl']], ['valid pathway'], [0.9])
    >>> validity_check.__globals__["logger"] = original_logger
    """

    valid_pathways: list[list[str]] = []
    valid_explanations: list[str] = []
    valid_confidence: list[float] = []

    for smiles_group, explanation, confidence in zip(
        res_molecules,
        res_explanations,
        res_confidence,
    ):
        candidate_pathway = (
            [smiles_group] if isinstance(smiles_group, str) else list(smiles_group)
        )
        valid_candidates: list[str] = []

        for smiles in candidate_pathway:
            if not is_valid_smiles(smiles):
                logger.warning(
                    "Invalid SMILES in pathway", molecule=molecule, smiles=smiles
                )
                continue

            if are_molecules_same(molecule, smiles):
                logger.warning(
                    "Molecule is same as target", molecule=molecule, smiles=smiles
                )
                continue

            if substructure_matching(molecule, smiles):
                logger.warning(
                    "Molecule is substructure of target",
                    molecule=molecule,
                    smiles=smiles,
                )
                continue

            valid_candidates.append(smiles)

        if len(valid_candidates) == len(candidate_pathway):
            valid_pathways.append(valid_candidates)
            valid_explanations.append(explanation)
            valid_confidence.append(confidence)

    logger.info("Validity check complete", valid_pathways=len(valid_pathways))
    return valid_pathways, valid_explanations, valid_confidence


def calc_mol_wt(mol: str) -> float:
    """Calculate the exact molecular weight for a SMILES string.

    Parameters
    ----------
    mol : str
        SMILES string of the molecule.

    Returns
    -------
    float
        Exact molecular weight. Returns ``0.0`` for invalid SMILES.

    Examples
    --------
    >>> round(calc_mol_wt("CCO"), 3)
    46.042
    >>> round(calc_mol_wt("C"), 3)
    16.031
    """

    molecule = _parse_molecule(mol)
    if molecule is None:
        logger.error("Error calculating molecular weight", smiles=mol)
        return 0.0

    return ExactMolWt(molecule)


def calc_chemical_formula(mol: str) -> str:
    """Calculate the molecular formula for a SMILES string.

    Parameters
    ----------
    mol : str
        SMILES string of the molecule.

    Returns
    -------
    str
        Molecular formula. Returns ``"N/A"`` for invalid SMILES.

    Examples
    --------
    >>> calc_chemical_formula("C")
    'CH4'
    >>> calc_chemical_formula("CCO")
    'C2H6O'
    """

    molecule = _parse_molecule(mol)
    if molecule is None:
        logger.error("Error calculating chemical formula", smiles=mol)
        return "N/A"

    return CalcMolFormula(molecule)


def are_molecules_same(smiles1: str, smiles2: str) -> bool:
    """Check whether two SMILES strings describe the same molecule.

    Parameters
    ----------
    smiles1 : str
        SMILES string of the first molecule.
    smiles2 : str
        SMILES string of the second molecule.

    Returns
    -------
    bool
        ``True`` when the molecules are equivalent, otherwise ``False``.

    Raises
    ------
    ValueError
        If either SMILES string is invalid.

    Examples
    --------
    >>> are_molecules_same("CCO", "OCC")
    True
    >>> are_molecules_same("CCO", "c1ccccc1")
    False
    """

    mol1 = _parse_molecule(smiles1)
    mol2 = _parse_molecule(smiles2)
    if mol1 is None or mol2 is None:
        raise ValueError(_INVALID_SMILES_ERROR)

    if Chem.MolToSmiles(mol1, canonical=True) == Chem.MolToSmiles(
        mol2,
        canonical=True,
    ):
        return True

    fingerprint1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol1,
        radius=2,
        nBits=1024,
    )
    fingerprint2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol2,
        radius=2,
        nBits=1024,
    )
    return fingerprint1 == fingerprint2


def compute_fingerprint(
    smiles: str,
    radius: int = 2,
    nBits: int = 2048,
) -> list[int] | None:
    """Compute a Morgan fingerprint for a molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    radius : int, optional
        Fingerprint radius, by default ``2``.
    nBits : int, optional
        Number of bits in the fingerprint, by default ``2048``.

    Returns
    -------
    list[int] | None
        Fingerprint bit vector as integers, or ``None`` when the SMILES is
        invalid.

    Examples
    --------
    >>> fingerprint = compute_fingerprint("CCO", radius=2, nBits=16)
    >>> len(fingerprint)
    16
    >>> compute_fingerprint("not_a_smiles") is None
    True
    """

    molecule = _parse_molecule(smiles)
    if molecule is None:
        return None

    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        molecule,
        radius,
        nBits=nBits,
    )
    return list(fingerprint)


def _has_ring_of_size(smiles: str, ring_size: int) -> bool:
    """Check whether a molecule contains a ring of the requested size."""

    molecule = _parse_molecule(smiles)
    if molecule is None:
        raise ValueError(_INVALID_SMILES_ERROR)

    return any(len(ring) == ring_size for ring in molecule.GetRingInfo().AtomRings())


def detect_seven_member_rings(smiles: str) -> bool:
    """Detect whether a molecule contains a seven-membered ring.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    bool
        ``True`` when at least one seven-membered ring is present.

    Raises
    ------
    ValueError
        If the SMILES string is invalid.

    Examples
    --------
    >>> detect_seven_member_rings("C1CCCCCC1")
    True
    >>> detect_seven_member_rings("C1CCCCC1")
    False
    """

    return _has_ring_of_size(smiles, ring_size=7)


def detect_eight_member_rings(smiles: str) -> bool:
    """Detect whether a molecule contains an eight-membered ring.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    bool
        ``True`` when at least one eight-membered ring is present.

    Raises
    ------
    ValueError
        If the SMILES string is invalid.

    Examples
    --------
    >>> detect_eight_member_rings("C1CCCCCCC1")
    True
    >>> detect_eight_member_rings("C1CCCCCC1")
    False
    """

    return _has_ring_of_size(smiles, ring_size=8)
