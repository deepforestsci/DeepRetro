"""Synthesis pathway visualization.

See :doc:`/package/deepretro.utils.visualize` for the module overview,
output layout, dependencies, and usage examples.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from io import BytesIO
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

TARGET_FILL = (210, 235, 210)
TARGET_STROKE = (130, 185, 130)
STEP_FILL = (210, 225, 245)
STEP_STROKE = (150, 180, 215)
BG = (255, 255, 255)
LINE_CLR = (180, 180, 180)
TEXT_CLR = (80, 80, 80)

MAIN_R = 105
SMALL_R = 68
COL_GAP = 310
MOL_CELL_H = 200
CHILD_GAP = 30
PAD = 80


@dataclass
class Node:
    """A single node in the rendered synthesis tree.

    Each node represents one synthesis step and holds the molecules
    displayed in that step's column plus layout coordinates assigned
    by :func:`layout`.

    Attributes
    ----------
    label : str
        Display label (e.g. ``"Step 1"``).
    molecules : list[dict]
        One dict per molecule rendered in the node, matching
        the ``reactants`` / ``products`` shape from
        :meth:`~deepretro.algorithms.autosolve.AutoSolver.solve`.
    is_target : bool
        ``True`` for the root node (drawn larger and green).
    children : list[Node]
        Child nodes (precursor steps) laid out to the right.
    x, y : int
        Pixel coordinates assigned during :func:`layout`.
    """

    label: str
    molecules: list[dict]
    is_target: bool = False
    children: list[Node] = field(default_factory=list)
    x: int = 0
    y: int = 0


def build_tree(result: dict[str, Any]) -> Node | None:
    """Build a :class:`Node` tree from an AutoSolver result.

    Walks the flat ``steps`` list plus ``dependencies`` map into a nested
    tree rooted at *Step 0* (the overall target molecule).

    Parameters
    ----------
    result : dict
        Output of :meth:`AutoSolver.solve`, containing ``"steps"`` and
        ``"dependencies"`` keys.

    Returns
    -------
    Node or None
        Root node of the tree, or ``None`` if *result* has no steps.

    Examples
    --------
    >>> from deepretro.utils.visualize import build_tree
    >>> result = {
    ...     "steps": [{"step": "1", "reactants": [{"smiles": "CCO"}],
    ...                "products": [{"smiles": "CC=O"}]}],
    ...     "dependencies": {"1": []},
    ... }
    >>> root = build_tree(result)
    >>> root.label
    'Step 0'
    """
    steps = result.get("steps", [])
    deps = result.get("dependencies", {})
    if not steps:
        return None

    step_map = {s["step"]: s for s in steps}
    step1 = step_map.get("1")
    if step1 is None:
        return None

    target = step1["products"][0] if step1.get("products") else None
    root = Node(
        label="Step 0",
        molecules=[target] if target else [],
        is_target=True,
    )

    def recurse(step_id: str) -> Node:
        """Build the subtree rooted at ``step_id``.

        Closes over ``step_map`` and ``deps`` from the enclosing
        :func:`build_tree` call.  The step's reactants become the
        node's molecules, and each dependency id listed in
        ``deps[step_id]`` is expanded into a child node via a
        recursive call.  Dependency ids pointing at steps absent
        from ``step_map`` are silently skipped so malformed inputs
        degrade gracefully rather than raising :class:`KeyError`.

        Parameters
        ----------
        step_id : str
            Key into ``step_map`` identifying the step to expand.

        Returns
        -------
        Node
            Fully populated subtree with label ``"Step <step_id>"``.
        """
        step = step_map[step_id]
        node = Node(
            label=f"Step {step_id}",
            molecules=step.get("reactants", []),
        )
        for child_id in deps.get(step_id, []):
            if child_id in step_map:
                node.children.append(recurse(child_id))
        return node

    root.children.append(recurse("1"))
    return root


def node_height(node: Node) -> int:
    """Compute total pixel height of a node's subtree.

    Used by :func:`layout` to evenly distribute vertical space between
    siblings without overlap.

    Parameters
    ----------
    node : Node
        Subtree root.

    Returns
    -------
    int
        Height in pixels.  At minimum ``MOL_CELL_H`` for a single
        molecule; larger when a node stacks multiple molecules or its
        subtree is taller than the node itself.

    Examples
    --------
    >>> from deepretro.utils.visualize import Node, node_height, MOL_CELL_H
    >>> leaf = Node(label="Step 1", molecules=[{"smiles": "CCO"}])
    >>> node_height(leaf) == MOL_CELL_H
    True
    """
    own = max(len(node.molecules), 1) * MOL_CELL_H
    if not node.children:
        return own
    kids = sum(node_height(c) for c in node.children) + \
        CHILD_GAP * (len(node.children) - 1)
    return max(own, kids)


def layout(node: Node, x: int, y_start: int, y_end: int) -> None:
    """Assign ``(x, y)`` coordinates to every node in the tree in-place.

    Parameters
    ----------
    node : Node
        Subtree root.  Mutated directly — the caller's tree now carries
        layout positions.
    x : int
        X coordinate for *node*.  Children are placed at ``x + COL_GAP``.
    y_start, y_end : int
        Vertical band allocated to *node*'s subtree.  The node is
        centered in this band; children are distributed within.

    Examples
    --------
    >>> from deepretro.utils.visualize import Node, layout
    >>> root = Node(label="Step 0", molecules=[], is_target=True)
    >>> layout(root, 100, 0, 400)
    >>> (root.x, root.y)
    (100, 200)
    """
    node.x = x
    node.y = (y_start + y_end) // 2

    if not node.children:
        return

    needed = sum(node_height(c) for c in node.children) + \
        CHILD_GAP * (len(node.children) - 1)
    start = y_start + max(0, (y_end - y_start - needed) // 2)

    for child in node.children:
        ch = node_height(child)
        layout(child, x + COL_GAP, start, start + ch)
        start += ch + CHILD_GAP


def render_mol(smiles: str, px: int) -> Image.Image | None:
    """Render a SMILES string into a transparent-background RGBA image.

    Falls back from the Cairo drawer to ``Draw.MolToImage`` if Cairo is
    unavailable.  White pixels are converted to transparent so the
    structure blends cleanly over the coloured node circles.

    Parameters
    ----------
    smiles : str
        Molecule to render.
    px : int
        Target image width/height in pixels (square).

    Returns
    -------
    PIL.Image.Image or None
        RGBA image, or ``None`` when RDKit cannot parse *smiles*.

    Examples
    --------
    >>> from deepretro.utils.visualize import render_mol  # doctest: +SKIP
    >>> img = render_mol("CCO", 200)                       # doctest: +SKIP
    >>> img.mode                                           # doctest: +SKIP
    'RGBA'
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(px, px)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.bondLineWidth = 2.0
        opts.padding = 0.15
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img = Image.open(BytesIO(drawer.GetDrawingText())).convert("RGBA")
    except Exception:
        img = Draw.MolToImage(mol, size=(px, px)).convert("RGBA")
    data = np.array(img)
    data[(data[:, :, :3] > 240).all(axis=2), 3] = 0
    return Image.fromarray(data)


def mol_metadata(mol_data: dict) -> tuple[str, str]:
    """Extract ``(formula, mass_string)`` from a molecule dict.

    Prefers pre-computed metadata (``product_metadata`` /
    ``reactant_metadata``) when present, otherwise recomputes via RDKit.

    Parameters
    ----------
    mol_data : dict
        Molecule entry from the AutoSolver output.  May contain a
        ``smiles`` key and/or a metadata dict.

    Returns
    -------
    tuple[str, str]
        ``(chemical_formula, "{mass} g/mol")``.  The mass string is
        empty when mass is zero/missing; ``"?"`` when RDKit cannot
        parse the SMILES.

    Examples
    --------
    >>> from deepretro.utils.visualize import mol_metadata  # doctest: +SKIP
    >>> mol_metadata({"smiles": "CCO"})                     # doctest: +SKIP
    ('C2H6O', '46.0 g/mol')
    >>> mol_metadata({"smiles": "INVALID"})                 # doctest: +SKIP
    ('INVALID', '?')
    """
    meta = mol_data.get("product_metadata") or mol_data.get("reactant_metadata")
    if meta:
        mass = meta.get("mass", 0)
        formula = meta.get("chemical_formula", "")
        return formula, f"{round(mass, 1)} g/mol" if mass else ""

    smiles = mol_data.get("smiles", "")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, "?"
    return rdMolDescriptors.CalcMolFormula(mol), f"{round(Descriptors.ExactMolWt(mol), 1)} g/mol"


def get_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Return a TrueType font, falling back to Pillow's default.

    Tries Arial then DejaVuSans (cross-platform).  If no TrueType file
    is available, returns Pillow's bitmap default so rendering still
    works in minimal environments.

    Parameters
    ----------
    size : int
        Desired font size in points.

    Returns
    -------
    PIL.ImageFont.FreeTypeFont or PIL.ImageFont.ImageFont

    Examples
    --------
    >>> from deepretro.utils.visualize import get_font  # doctest: +SKIP
    >>> font = get_font(14)                              # doctest: +SKIP
    """
    for name in ("arial.ttf", "Arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def text_centered(
    draw: ImageDraw.ImageDraw,
    cx: int,
    y: int,
    text: str,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    fill: tuple[int, int, int] = TEXT_CLR,
) -> None:
    """Draw *text* horizontally centered at ``(cx, y)``.

    Parameters
    ----------
    draw : PIL.ImageDraw.ImageDraw
        Active draw context.
    cx : int
        Target center x.  The text's left edge is offset so its midpoint
        sits on *cx*.
    y : int
        Top y coordinate for the text.
    text : str
        String to render.
    font : PIL.ImageFont.FreeTypeFont or PIL.ImageFont.ImageFont
        Font used for measurement and rendering.
    fill : tuple[int, int, int]
        RGB colour, default :data:`TEXT_CLR`.
    """
    bb = draw.textbbox((0, 0), text, font=font)
    draw.text((cx - (bb[2] - bb[0]) // 2, y), text, fill=fill, font=font)


def max_x(node: Node) -> int:
    """Return the rightmost x coordinate anywhere in the tree.

    Used to size the output canvas wide enough to contain every child.

    Parameters
    ----------
    node : Node
        Subtree root (typically the overall root after :func:`layout`).

    Returns
    -------
    int
        Largest ``x`` attribute found.

    Examples
    --------
    >>> from deepretro.utils.visualize import Node, max_x
    >>> root = Node(label="Step 0", molecules=[])
    >>> root.x = 100
    >>> max_x(root)
    100
    """
    result = node.x
    for c in node.children:
        result = max(result, max_x(c))
    return result


def draw_edges(draw: ImageDraw.ImageDraw, node: Node) -> None:
    """Draw cubic bezier connectors from *node* to each descendant.

    Recursively traverses the tree so every parent-child edge in the
    subtree rooted at *node* is drawn once.

    Parameters
    ----------
    draw : PIL.ImageDraw.ImageDraw
        Active draw context on the output canvas.
    node : Node
        Subtree root.  Must have layout coordinates already assigned.
    """
    for child in node.children:
        sx, sy = node.x, node.y
        tx, ty = child.x, child.y
        mid = (sx + tx) // 2
        pts: list[tuple[int, int]] = []
        for i in range(31):
            t = i / 30
            u = 1 - t
            bx = u**3 * sx + 3 * u**2 * t * mid + 3 * u * t**2 * mid + t**3 * tx
            by = u**3 * sy + 3 * u**2 * t * sy + 3 * u * t**2 * ty + t**3 * ty
            pts.append((int(bx), int(by)))
        for j in range(len(pts) - 1):
            draw.line([pts[j], pts[j + 1]], fill=LINE_CLR, width=2)
        draw_edges(draw, child)


def draw_node(
    img: Image.Image,
    draw: ImageDraw.ImageDraw,
    node: Node,
    font_lbl: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    font_sm: ImageFont.FreeTypeFont | ImageFont.ImageFont,
) -> None:
    """Render step label, molecule circles, structures, and metadata.

    Walks the tree depth-first.  For each node:

    1. Draw the step label above the column.
    2. Draw one circle per molecule (green if target, blue otherwise).
    3. Paste the rendered structure inside each circle.
    4. Add formula and mass below each circle.

    Parameters
    ----------
    img : PIL.Image.Image
        Canvas to paste molecule structures onto.
    draw : PIL.ImageDraw.ImageDraw
        Draw context for the same canvas.
    node : Node
        Subtree root (with layout already applied).
    font_lbl : PIL.ImageFont
        Font for step labels.
    font_sm : PIL.ImageFont
        Font for formula and mass annotations.
    """
    mols = node.molecules
    n = max(len(mols), 1)
    r = MAIN_R if node.is_target else SMALL_R
    top_y = node.y - (n * MOL_CELL_H) // 2

    text_centered(draw, node.x, top_y - 24, node.label, font_lbl)

    for i, mol_data in enumerate(mols):
        cy = top_y + i * MOL_CELL_H + MOL_CELL_H // 2
        cx = node.x

        fill = TARGET_FILL if node.is_target else STEP_FILL
        stroke = TARGET_STROKE if node.is_target else STEP_STROKE
        draw.ellipse(
            [cx - r, cy - r, cx + r, cy + r],
            fill=fill, outline=stroke, width=2,
        )

        smiles = mol_data.get("smiles", "") if isinstance(mol_data, dict) else str(mol_data)
        mol_px = int(r * 1.5)
        mol_img = render_mol(smiles, mol_px + 40)
        if mol_img is not None:
            mol_img = mol_img.resize((mol_px, mol_px), Image.Resampling.LANCZOS)
            img.paste(mol_img, (cx - mol_px // 2, cy - mol_px // 2), mol_img)

        if isinstance(mol_data, dict):
            formula, mass_str = mol_metadata(mol_data)
        else:
            formula, mass_str = str(mol_data), ""

        text_centered(draw, cx, cy + r + 6, mass_str, font_sm)
        text_centered(draw, cx, cy + r + 22, formula, font_sm)

    for child in node.children:
        draw_node(img, draw, child, font_lbl, font_sm)


def visualize_pathway(result: dict[str, Any]) -> Image.Image:
    """Render a retrosynthesis pathway as a PIL image.

    Top-level entry point.  Builds the tree, lays it out, sizes the
    canvas, draws bezier edges, then draws every node.

    Parameters
    ----------
    result : dict
        Output of :meth:`AutoSolver.solve`, with ``"steps"`` and
        ``"dependencies"`` keys.

    Returns
    -------
    PIL.Image.Image
        RGB image of the full pathway.  When *result* is empty a small
        placeholder image saying "No synthesis steps" is returned.

    Examples
    --------
    >>> from deepretro.algorithms.autosolve import AutoSolver  # doctest: +SKIP
    >>> from deepretro.utils.visualize import visualize_pathway  # doctest: +SKIP
    >>> solver = AutoSolver()                                  # doctest: +SKIP
    >>> result = solver.solve("CCO")                           # doctest: +SKIP
    >>> img = visualize_pathway(result)                        # doctest: +SKIP
    >>> img.save("out.png")                                    # doctest: +SKIP
    """
    font_lbl = get_font(15)
    font_sm = get_font(13)

    root = build_tree(result)
    if root is None:
        img = Image.new("RGB", (400, 200), BG)
        draw = ImageDraw.Draw(img)
        text_centered(draw, 200, 90, "No synthesis steps", font_lbl)
        return img

    tree_h = node_height(root)
    canvas_h = tree_h + PAD * 2
    layout(root, PAD + MAIN_R, PAD, PAD + tree_h)

    canvas_w = max_x(root) + MAIN_R + PAD + 60

    img = Image.new("RGBA", (canvas_w, canvas_h), (*BG, 255))
    draw = ImageDraw.Draw(img)

    draw_edges(draw, root)
    draw_node(img, draw, root, font_lbl, font_sm)

    return img.convert("RGB")
