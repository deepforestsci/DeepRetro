"""Synthesis pathway visualization.

Mirrors the viewer's (app_v4.js) processData + renderGraph logic as a
static PIL Image:

- Step 0 = target molecule (product of step 1) — green circle
- Step N = reactants of that step — blue circles, stacked vertically
- Layout follows the full dependency tree left-to-right.

Dependencies: PIL (Pillow), numpy, rdkit.
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
    label: str
    molecules: list[dict]
    is_target: bool = False
    children: list[Node] = field(default_factory=list)
    x: int = 0
    y: int = 0


def build_tree(result: dict[str, Any]) -> Node | None:
    """Build a tree of Nodes from autosolve output, mirroring the JS viewer."""
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
    """Pixel height needed for this node's subtree."""
    own = max(len(node.molecules), 1) * MOL_CELL_H
    if not node.children:
        return own
    kids = sum(node_height(c) for c in node.children) + \
           CHILD_GAP * (len(node.children) - 1)
    return max(own, kids)


def layout(node: Node, x: int, y_start: int, y_end: int):
    """Assign (x, y) to every node in the tree."""
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
    """SMILES → transparent-background RGBA image."""
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
    """Extract (formula, mass_str) from a molecule dict."""
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


def get_font(size: int):
    for name in ("arial.ttf", "Arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def text_centered(draw, cx, y, text, font, fill=TEXT_CLR):
    bb = draw.textbbox((0, 0), text, font=font)
    draw.text((cx - (bb[2] - bb[0]) // 2, y), text, fill=fill, font=font)


def max_x(node: Node) -> int:
    """Rightmost x coordinate in the tree."""
    result = node.x
    for c in node.children:
        result = max(result, max_x(c))
    return result


def draw_edges(draw: ImageDraw.ImageDraw, node: Node):
    """Draw cubic bezier connectors from parent to each child."""
    for child in node.children:
        sx, sy = node.x, node.y
        tx, ty = child.x, child.y
        mid = (sx + tx) // 2
        pts = []
        for i in range(31):
            t = i / 30
            u = 1 - t
            bx = u**3*sx + 3*u**2*t*mid + 3*u*t**2*mid + t**3*tx
            by = u**3*sy + 3*u**2*t*sy  + 3*u*t**2*ty  + t**3*ty
            pts.append((int(bx), int(by)))
        for j in range(len(pts) - 1):
            draw.line([pts[j], pts[j + 1]], fill=LINE_CLR, width=2)
        draw_edges(draw, child)


def draw_node(img: Image.Image, draw: ImageDraw.ImageDraw, node: Node,
              font_lbl, font_sm):
    """Draw step label, molecule circles, structures, and metadata."""
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
        draw.ellipse([cx - r, cy - r, cx + r, cy + r],
                     fill=fill, outline=stroke, width=2)

        smiles = mol_data.get("smiles", "") if isinstance(mol_data, dict) else str(mol_data)
        mol_px = int(r * 1.5)
        mol_img = render_mol(smiles, mol_px + 40)
        if mol_img is not None:
            mol_img = mol_img.resize((mol_px, mol_px), Image.LANCZOS)
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
    """Render a retrosynthesis pathway as a PIL Image.

    Parameters
    ----------
    result : dict
        Output of ``autosolve()``, with ``"steps"`` and ``"dependencies"``.

    Returns
    -------
    PIL.Image.Image
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
