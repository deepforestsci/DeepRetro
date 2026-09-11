"""Reversible abstraction and chemical deprotection structure candidates.

Motif matches are candidates, not evidence of a group's synthetic purpose.
The retained O/N atom remains part of the scaffold; only the pendant group
is hidden in masking mode. Proposal mode generates one-site deprotected product
structures; it does not establish reaction conditions or chemical feasibility.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
from uuid import uuid4

from rdkit import Chem

# Larger motifs precede their contained smaller motifs. Atom 0 is the retained
# scaffold heteroatom; atom 1 is the attachment atom of the hidden fragment.
_PATTERNS = {
    "Boc": "[N;X3]-C(=O)-O-C(C)(C)C",
    "Cbz": "[N;X3]-C(=O)-O-[CH2]-c1ccccc1",
    "TBS": "[O;X2]-[Si]([CH3])([CH3])-C([CH3])([CH3])[CH3]",
    "TES": "[O;X2]-[Si]([CH2][CH3])([CH2][CH3])[CH2][CH3]",
    "TBDPS": "[O;X2]-[Si](c1ccccc1)(c1ccccc1)C([CH3])([CH3])[CH3]",
    "OBn": "[O;X2]-[CH2]-c1ccccc1",
    "OEt": "[O;X2]-[CH2]-[CH3]",
    "OMe": "[O;X2]-[CH3]",
}
SUPPORTED_GROUPS = tuple(_PATTERNS)

# Every pattern above anchors on an O whose *other* substituent is left
# unconstrained by the SMARTS itself — intentionally, so the same pattern
# matches the anchor sitting on any carbon scaffold. But that also means an
# anchor O bonded to something else entirely (e.g. a phosphonate's P-O-CH3)
# satisfies the SMARTS just as well as a real alkyl ether does, even though
# it isn't a protecting group at all — it's a reactive handle (see
# architecture discussion: HWE phosphonate esters). Enforced separately below
# rather than folded into each SMARTS, since it's the same rule for all of
# them: the scaffold-side neighbor of an oxygen anchor must be carbon.
_OXYGEN_ANCHOR_GROUPS = {"TBS", "TES", "TBDPS", "OBn", "OEt", "OMe"}


@dataclass(frozen=True)
class _ProtectionSite:
    """One accepted pendant motif in a canonical molecule's atom ordering."""

    group: str
    anchor: int
    hidden: tuple[int, ...]


def _parse(smiles: str) -> Any:
    """Parse nonempty molecular SMILES without silently accepting names."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a nonempty string")
    params = Chem.SmilesParserParams()
    params.parseName = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("invalid SMILES")
    return mol


def _find_sites(
    smiles: str, groups: list[str] | None
) -> tuple[Any, list[_ProtectionSite]]:
    """Canonicalize full input and identify supported, nonoverlapping sites."""
    mol = _parse(smiles)
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ValueError("provide a full molecule without dummy atoms")
    mol = _parse(Chem.MolToSmiles(mol, isomericSmiles=True))
    if groups is not None and (
        not isinstance(groups, list)
        or any(not isinstance(g, str) or g not in _PATTERNS for g in groups)
    ):
        raise ValueError(f"groups must be a list drawn from {list(SUPPORTED_GROUPS)}")
    selected = set(SUPPORTED_GROUPS if groups is None else groups)
    used: set[int] = set()
    sites: list[_ProtectionSite] = []
    for name, smarts in _PATTERNS.items():
        if name not in selected:
            continue
        for match in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
            anchor = match[0]
            hidden = set(match[1:])
            if used.intersection(match):
                continue
            # Rebuilding a chiral fragment changes neighbor order. Leave it
            # explicit rather than risk stereo inversion during restoration.
            if any(
                mol.GetAtomWithIdx(idx).GetChiralTag()
                != Chem.ChiralType.CHI_UNSPECIFIED
                for idx in hidden
            ):
                continue
            # Exclude terminal OH and fragments with multiple attachments.
            anchor_atom = mol.GetAtomWithIdx(anchor)
            if anchor_atom.GetDegree() < 2:
                continue
            # An oxygen-anchored pattern (ether/silyl ether) is only a real
            # protecting group when the scaffold side of the anchor is
            # carbon. A P-O-CH3 (phosphonate/phosphate ester) or S-O-CH3
            # (sulfonate ester) satisfies the same SMARTS but is a reactive
            # functional group, not a masked alcohol — see _OXYGEN_ANCHOR_GROUPS.
            if name in _OXYGEN_ANCHOR_GROUPS and any(
                neighbor.GetIdx() not in hidden and neighbor.GetSymbol() != "C"
                for neighbor in anchor_atom.GetNeighbors()
            ):
                continue
            boundary = [
                bond
                for bond in mol.GetBonds()
                if (bond.GetBeginAtomIdx() in hidden)
                != (bond.GetEndAtomIdx() in hidden)
            ]
            if len(boundary) != 1 or boundary[0].GetBondType() != Chem.BondType.SINGLE:
                continue
            sites.append(_ProtectionSite(name, anchor, tuple(match[1:])))
            used.update(match)
    return mol, sites


class ProtectionContext:
    """Keep exact hidden fragments private to one agent tool registry.

    Examples
    --------
    >>> context = ProtectionContext()
    >>> result = context.handle_protection("COc1ccccc1", groups=["OMe"])
    >>> context.handle_deprotection(result["masked_smiles"], result["mask_id"])["smiles"]
    'COc1ccccc1'
    """

    def __init__(self) -> None:
        self._masks: dict[str, dict[int, tuple[Any, _ProtectionSite]]] = {}

    def handle_protection(
        self, smiles: str, groups: list[str] | None = None
    ) -> dict[str, Any]:
        """Hide recognized pendant groups behind mapped dummy atoms.

        Parameters
        ----------
        smiles : str
            Full molecular SMILES, including any salts and atom maps.
        groups : list of str or None, optional
            Motifs to consider: Boc, Cbz, TBS, OBn, OEt, OMe. Defaults to all.
            An empty list disables masking.

        Returns
        -------
        dict
            Canonical original, masked SMILES, mask ID, matched group legend,
            and reasoning guidance. Restoration data stays in this context.
        """
        mol, sites = _find_sites(smiles, groups)
        original = Chem.MolToSmiles(mol, isomericSmiles=True)
        editable = Chem.RWMol(mol)
        remove: set[int] = set()
        fragments: dict[int, tuple[Any, _ProtectionSite]] = {}
        legend: list[dict[str, Any]] = []
        first_label = max(atom.GetAtomMapNum() for atom in mol.GetAtoms()) + 1
        for label, site in enumerate(sites, start=first_label):
            root, *rest = site.hidden
            dummy = Chem.Atom(0)
            dummy.SetAtomMapNum(label)
            editable.ReplaceAtom(root, dummy)
            fragments[label] = (mol, site)
            remove.update(rest)
            legend.append({"group": site.group, "marker": f"[*:{label}]"})
        for idx in sorted(remove, reverse=True):
            editable.RemoveAtom(idx)
        masked = editable.GetMol()
        Chem.SanitizeMol(masked)
        mask_id = f"mask_{uuid4().hex}"
        self._masks[mask_id] = fragments
        return {
            "original_smiles": original,
            "masked_smiles": Chem.MolToSmiles(masked, isomericSmiles=True),
            "mask_id": mask_id,
            "groups": legend,
            "matched": bool(legend),
            "guidance": (
                "These are candidate protecting-group motifs, not confirmed roles. "
                "Focus retrosynthetic disconnections on the exposed core. Keep each "
                "mapped dummy and its attachment unchanged; consider compatibility "
                "with the hidden groups. Call handle_deprotection on each masked "
                "precursor with this mask_id before validation and final output. "
                "Use only full, restored molecular SMILES in the final answer."
            ),
        }

    def handle_deprotection(
        self,
        smiles: str,
        mask_id: str | None = None,
        mode: str = "restore",
        groups: list[str] | None = None,
    ) -> dict[str, Any]:
        """Restore masked groups or propose chemical deprotection products.

        Parameters
        ----------
        smiles : str
            Masked SMILES for restoration; full protected SMILES for proposals.
        mask_id : str or None, optional
            Required for restoration: ID from handle_protection in this run.
            Must be omitted for proposals.
        mode : {"restore", "propose"}, optional
            Restore hidden atoms (default), or propose one-site forward chemical
            deprotection products without altering the restoration context.
        groups : list of str or None, optional
            Proposal-only filter over the supported motifs; None selects all.

        Returns
        -------
        dict
            Restoration returns full ``smiles`` and ``restored_groups``.
            Proposal returns full product ``candidates`` with explicit forward
            direction and unverified reaction conditions.

        Examples
        --------
        >>> context = ProtectionContext()
        >>> result = context.handle_deprotection(
        ...     "CCNC(=O)OC(C)(C)C", mode="propose", groups=["Boc"]
        ... )
        >>> result["candidates"][0]["product_smiles"]
        'CCN'
        """
        if mode == "propose":
            if mask_id is not None:
                raise ValueError("propose mode requires full SMILES and no mask_id")
            return self._propose_deprotection(smiles, groups)
        if mode != "restore":
            raise ValueError("mode must be 'restore' or 'propose'")
        if groups is not None:
            raise ValueError("groups is only supported in propose mode")
        if not isinstance(mask_id, str) or mask_id not in self._masks:
            raise ValueError("unknown mask_id; use one from handle_protection")
        mol = _parse(smiles)
        editable = Chem.RWMol(mol)
        fragments = self._masks[mask_id]
        seen: set[int] = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 0:
                continue
            label = atom.GetAtomMapNum()
            if label not in fragments or label in seen:
                raise ValueError("unknown or duplicate protecting-group marker")
            if atom.GetDegree() != 1 or atom.GetIsotope() or atom.GetFormalCharge():
                raise ValueError("a marker must have one unchanged single attachment")
            bond = atom.GetBonds()[0]
            source, site = fragments[label]
            indices = site.hidden
            root = indices[0]
            original_anchor = source.GetAtomWithIdx(site.anchor)
            if (
                bond.GetBondType() != Chem.BondType.SINGLE
                or atom.GetNeighbors()[0].GetAtomicNum()
                != original_anchor.GetAtomicNum()
            ):
                raise ValueError("marker attachment must retain its original O/N atom")
            idx = atom.GetIdx()
            editable.ReplaceAtom(idx, Chem.Atom(source.GetAtomWithIdx(root)))
            mapping = {root: idx}
            for old in indices[1:]:
                mapping[old] = editable.AddAtom(Chem.Atom(source.GetAtomWithIdx(old)))
            for original_bond in source.GetBonds():
                begin, end = (
                    original_bond.GetBeginAtomIdx(),
                    original_bond.GetEndAtomIdx(),
                )
                if begin in mapping and end in mapping:
                    editable.AddBond(
                        mapping[begin], mapping[end], original_bond.GetBondType()
                    )
            seen.add(label)
        restored = editable.GetMol()
        Chem.SanitizeMol(restored)
        return {
            "mode": "restore",
            "smiles": Chem.MolToSmiles(restored, isomericSmiles=True),
            "restored_groups": len(seen),
            "guidance": "Groups restored; no chemical deprotection reaction was performed.",
        }

    def _propose_deprotection(
        self, smiles: str, groups: list[str] | None
    ) -> dict[str, Any]:
        """Generate one-site cleavage candidates using the same motif matcher."""
        source, sites = _find_sites(smiles, groups)
        original = Chem.MolToSmiles(source, isomericSmiles=True)
        candidates: list[dict[str, Any]] = []
        for site in sites:
            hidden = site.hidden
            anchor = source.GetAtomWithIdx(site.anchor)
            # Cleavage changes this atom's substituent order. Conservatively
            # exclude unusual stereogenic O/N anchors rather than invert them.
            if anchor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            editable = Chem.RWMol(source)
            retained = editable.GetAtomWithIdx(anchor.GetIdx())
            retained.SetNumExplicitHs(anchor.GetTotalNumHs() + 1)
            retained.SetNoImplicit(True)
            for idx in sorted(hidden, reverse=True):
                editable.RemoveAtom(idx)
            product = editable.GetMol()
            Chem.SanitizeMol(product)
            product_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
            candidates.append(
                {
                    "protected_smiles": original,
                    "product_smiles": product_smiles,
                    "removed_group": site.group,
                    "site_atom_index": anchor.GetIdx(),
                    "site_atom_map": anchor.GetAtomMapNum(),
                    "reaction_smiles": f"{original}>>{product_smiles}",
                    "conditions": None,
                    "status": "structural_candidate",
                }
            )
        return {
            "mode": "propose",
            "direction": "forward",
            "candidates": candidates,
            "guidance": (
                "Each candidate removes one matched motif and restores an O-H or "
                "N-H bond; other sites remain protected. These are forward "
                "protected-substrate -> deprotected-product proposals, not "
                "retrosynthetic precursors of the protected input. A motif may "
                "be an intended substituent rather than a protecting group. "
                "Conditions, selectivity, and substrate compatibility are not "
                "validated; use the full structures to assess them before "
                "proposing a chemical step. Reaction SMILES omit reagents and "
                "byproducts and are not atom-balanced. An empty list means no "
                "supported site was found, not that deprotection is impossible."
            ),
        }
