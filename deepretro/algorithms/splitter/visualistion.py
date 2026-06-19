"""
Visualization utilities for r-BRICS molecule splitting.

Examples
--------
>>> from deepretro.algorithms.splitter.rbrics_splitter import split_molecule
>>> from deepretro.algorithms.splitter.rbrics_viz import draw_split
>>> fragments = split_molecule(smiles)
>>> img = draw_split(smiles, fragments)
>>> img.save("split.png")
"""

from typing import List, Optional, Tuple
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

from .rBRICS_public import FindrBRICSBonds


# Fragment colors (RGBA)
FRAG_COLORS = [
    (0.22, 0.56, 0.87, 0.35),  # blue
    (0.90, 0.30, 0.24, 0.35),  # red
    (0.30, 0.75, 0.40, 0.35),  # green
    (0.80, 0.55, 0.10, 0.35),  # amber
    (0.58, 0.30, 0.75, 0.35),  # purple
    (0.10, 0.70, 0.70, 0.35),  # teal
    (0.85, 0.45, 0.65, 0.35),  # pink
]


def draw_split(
    smiles: str,
    fragments: List[str],
    img_size: Tuple[int, int] = (400, 300),
    highlight_bonds: bool = True
) -> "PIL.Image.Image":
    """
    Draw the original molecule with fragments highlighted.

    Parameters
    ----------
    smiles : str
        Original molecule SMILES.
    fragments : List[str]
        Fragment SMILES from split_molecule().
    img_size : Tuple[int, int], default (400, 300)
        Image dimensions (width, height) in pixels.
    highlight_bonds : bool, default True
        If True, highlight the bonds that were cut in red.

    Returns
    -------
    PIL.Image.Image
        Image of the molecule with cut bonds highlighted.

    Examples
    --------
    >>> img = draw_split("CCCOCc1ccccc1", ["[4*]CCC", "[3*]O[3*]", "[4*]Cc1ccccc1"])
    >>> img.save("molecule_split.png")
    """
    from PIL import Image

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    AllChem.Compute2DCoords(mol)

    highlight_atoms = []
    highlight_bonds_list = []
    atom_colors = {}
    bond_colors = {}

    if highlight_bonds:
        cut_bonds = _find_cut_bonds(mol)
        highlight_bonds_list = cut_bonds
        bond_colors = {b: (0.9, 0.15, 0.15) for b in cut_bonds}

        for bid in cut_bonds:
            bond = mol.GetBondWithIdx(bid)
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            highlight_atoms.extend([a1, a2])
            atom_colors[a1] = (0.9, 0.15, 0.15)
            atom_colors[a2] = (0.9, 0.15, 0.15)

    d = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    opts = d.drawOptions()
    opts.useBWAtomPalette()
    opts.bondLineWidth = 2.5

    d.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=highlight_bonds_list,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    d.FinishDrawing()

    return Image.open(BytesIO(d.GetDrawingText()))


def draw_fragments(
    fragments: List[str],
    img_size: Tuple[int, int] = (250, 200),
    mols_per_row: int = 4
) -> "PIL.Image.Image":
    """
    Draw all fragments in a grid.

    Parameters
    ----------
    fragments : List[str]
        Fragment SMILES from split_molecule().
    img_size : Tuple[int, int], default (250, 200)
        Size of each fragment image.
    mols_per_row : int, default 4
        Number of molecules per row in the grid.

    Returns
    -------
    PIL.Image.Image
        Grid image of all fragments.

    Examples
    --------
    >>> img = draw_fragments(["[4*]CCC", "[3*]O[3*]", "[4*]Cc1ccccc1"])
    >>> img.save("fragments.png")
    """
    mols = []
    legends = []

    for i, smi in enumerate(fragments):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
            n_atoms = mol.GetNumHeavyAtoms()
            legends.append(f"Frag {i} ({n_atoms} atoms)")

    if not mols:
        raise ValueError("No valid fragments to draw")

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=img_size,
        legends=legends,
    )

    return img


def draw_split_summary(
    smiles: str,
    fragments: List[str],
    save_path: Optional[str] = None
) -> "PIL.Image.Image":
    """
    Draw complete split summary: original molecule + fragments grid.

    Parameters
    ----------
    smiles : str
        Original molecule SMILES.
    fragments : List[str]
        Fragment SMILES from split_molecule().
    save_path : str, optional
        If provided, save the image to this path.

    Returns
    -------
    PIL.Image.Image
        Combined image with original molecule on top and fragments below.

    Examples
    --------
    >>> img = draw_split_summary(smiles, fragments, save_path="summary.png")
    """
    from PIL import Image

    mol_img = draw_split(smiles, fragments, img_size=(500, 350))

    if len(fragments) > 1:
        frag_img = draw_fragments(fragments, img_size=(200, 150))
    else:
        frag_img = None

    if frag_img is None:
        combined = mol_img
    else:
        total_height = mol_img.height + frag_img.height + 10
        total_width = max(mol_img.width, frag_img.width)

        combined = Image.new("RGB", (total_width, total_height), "white")
        combined.paste(mol_img, ((total_width - mol_img.width) // 2, 0))
        combined.paste(frag_img, ((total_width - frag_img.width) // 2, mol_img.height + 10))

    if save_path:
        combined.save(save_path)

    return combined


def draw_full_pipeline(
    smiles: str,
    fragments: List[str],
    reactive_fragments: List[str],
    info: List[dict],
    save_path: Optional[str] = None,
) -> "PIL.Image.Image":
    """Draw the complete split-and-cap pipeline as a single image.

    Layout (top to bottom):
    1. Original molecule with cut bonds highlighted
    2. Raw fragments (with ``[n*]`` dummies)
    3. Reactive fragments (dummies replaced with caps)

    Parameters
    ----------
    smiles : str
        Original molecule SMILES.
    fragments : List[str]
        Raw fragment SMILES from ``split_molecule``.
    reactive_fragments : List[str]
        Capped fragment SMILES from ``replace_dummy_atoms``.
    info : List[dict]
        Replacement log from ``replace_dummy_atoms``.
    save_path : str, optional
        If provided, save the image to this path.

    Returns
    -------
    PIL.Image.Image
    """
    from PIL import Image, ImageDraw, ImageFont

    mol_img = draw_split(smiles, fragments, img_size=(600, 350))
    raw_img = _draw_labeled_grid(fragments, "Raw fragments (dummy atoms)")
    cap_img = _draw_labeled_grid(reactive_fragments, "Reactive fragments (capped)")

    padding = 10
    total_width = max(mol_img.width, raw_img.width, cap_img.width) + 2 * padding
    total_height = mol_img.height + raw_img.height + cap_img.height + 4 * padding

    canvas = Image.new("RGB", (total_width, total_height), "white")
    y = padding
    canvas.paste(mol_img, ((total_width - mol_img.width) // 2, y))
    y += mol_img.height + padding
    canvas.paste(raw_img, ((total_width - raw_img.width) // 2, y))
    y += raw_img.height + padding
    canvas.paste(cap_img, ((total_width - cap_img.width) // 2, y))

    if save_path:
        canvas.save(save_path)

    return canvas


def _draw_labeled_grid(
    fragments: List[str],
    title: str,
    img_size: Tuple[int, int] = (250, 200),
    mols_per_row: int = 4,
) -> "PIL.Image.Image":
    """Draw a fragment grid with a title banner."""
    from PIL import Image, ImageDraw

    mols = []
    legends = []
    for i, smi in enumerate(fragments):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
            legends.append(smi if len(smi) < 30 else f"Frag {i}")

    if not mols:
        img = Image.new("RGB", (400, 50), "white")
        d = ImageDraw.Draw(img)
        d.text((10, 15), "No valid fragments", fill="gray")
        return img

    grid = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=img_size,
        legends=legends,
    )
    if not isinstance(grid, Image.Image):
        from io import BytesIO
        grid = Image.open(BytesIO(grid.data))

    banner_h = 30
    combined = Image.new("RGB", (grid.width, grid.height + banner_h), "white")
    draw = ImageDraw.Draw(combined)
    draw.text((10, 5), title, fill="black")
    combined.paste(grid, (0, banner_h))
    return combined


def _find_cut_bonds(mol: Chem.Mol) -> List[int]:
    """Find bond indices that r-BRICS would cut."""
    cut_bonds = []
    seen = set()

    for (a1, a2), (t1, t2) in FindrBRICSBonds(mol):
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is not None:
            bid = bond.GetIdx()
            if bid not in seen:
                seen.add(bid)
                cut_bonds.append(bid)

    return cut_bonds