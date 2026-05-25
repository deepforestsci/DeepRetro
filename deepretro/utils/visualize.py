"""Render formatted retrosynthesis pathways as static images."""

from __future__ import annotations

from dataclasses import dataclass, field
from io import BytesIO
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

BACKGROUND = (255, 255, 255)
LINE_COLOR = (180, 180, 180)
TEXT_COLOR = (80, 80, 80)
TARGET_FILL = (210, 235, 210)
TARGET_STROKE = (130, 185, 130)
STEP_FILL = (210, 225, 245)
STEP_STROKE = (150, 180, 215)
MAIN_RADIUS = 105
SMALL_RADIUS = 68
COLUMN_GAP = 310
MOLECULE_CELL_HEIGHT = 200
CHILD_GAP = 30
PADDING = 80


@dataclass
class PathwayNode:
    """Renderable pathway node."""

    label: str
    molecules: list[dict[str, Any]]
    is_target: bool = False
    children: list["PathwayNode"] = field(default_factory=list)
    x: int = 0
    y: int = 0


def build_tree(result: dict[str, Any]) -> PathwayNode | None:
    """Build a render tree from formatted AutoSolver output."""
    steps = result.get("steps", [])
    dependencies = result.get("dependencies", {})
    if not steps:
        return None

    step_map = {step["step"]: step for step in steps}
    first_step = step_map.get("1")
    if first_step is None:
        return None

    target = first_step["products"][0] if first_step.get("products") else None
    root = PathwayNode(
        label="Step 0",
        molecules=[target] if target else [],
        is_target=True,
    )

    def _recurse(step_id: str) -> PathwayNode:
        step = step_map[step_id]
        node = PathwayNode(
            label=f"Step {step_id}",
            molecules=step.get("reactants", []),
        )
        for child_id in dependencies.get(step_id, []):
            if child_id in step_map:
                node.children.append(_recurse(child_id))
        return node

    root.children.append(_recurse("1"))
    return root


def node_height(node: PathwayNode) -> int:
    """Return the pixel height needed for a node subtree."""
    own_height = max(len(node.molecules), 1) * MOLECULE_CELL_HEIGHT
    if not node.children:
        return own_height
    child_height = sum(node_height(child) for child in node.children)
    child_height += CHILD_GAP * (len(node.children) - 1)
    return max(own_height, child_height)


def layout(node: PathwayNode, x: int, y_start: int, y_end: int) -> None:
    """Assign coordinates to every node in a render tree."""
    node.x = x
    node.y = (y_start + y_end) // 2
    if not node.children:
        return

    needed_height = sum(node_height(child) for child in node.children)
    needed_height += CHILD_GAP * (len(node.children) - 1)
    y = y_start + max(0, (y_end - y_start - needed_height) // 2)
    for child in node.children:
        child_height = node_height(child)
        layout(child, x + COLUMN_GAP, y, y + child_height)
        y += child_height + CHILD_GAP


def render_molecule(smiles: str, size: int) -> Image.Image | None:
    """Render one SMILES string to a transparent molecule image."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None

    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
        options = drawer.drawOptions()
        options.clearBackground = True
        options.bondLineWidth = 2.0
        options.padding = 0.15
        drawer.DrawMolecule(molecule)
        drawer.FinishDrawing()
        image = Image.open(BytesIO(drawer.GetDrawingText())).convert("RGBA")
    except Exception:
        image = Draw.MolToImage(molecule, size=(size, size)).convert("RGBA")

    data = np.array(image)
    data[(data[:, :, :3] > 240).all(axis=2), 3] = 0
    return Image.fromarray(data)


def molecule_metadata(molecule_data: dict[str, Any]) -> tuple[str, str]:
    """Return formula and formatted mass for a molecule entry."""
    metadata = molecule_data.get("product_metadata") or molecule_data.get(
        "reactant_metadata"
    )
    if metadata:
        mass = metadata.get("mass", 0)
        formula = str(metadata.get("chemical_formula", ""))
        return formula, f"{round(float(mass), 1)} g/mol" if mass else ""

    smiles = str(molecule_data.get("smiles", ""))
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return smiles, "?"
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    mass = round(Descriptors.ExactMolWt(molecule), 1)
    return formula, f"{mass} g/mol"


def get_font(size: int) -> Any:
    """Load a common TrueType font or fall back to PIL's default font."""
    for font_name in ("Arial.ttf", "arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(font_name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def draw_centered_text(
    draw: ImageDraw.ImageDraw,
    center_x: int,
    y: int,
    text: str,
    font: ImageFont.ImageFont,
) -> None:
    """Draw centered text at a given y coordinate."""
    bbox = draw.textbbox((0, 0), text, font=font)
    width = bbox[2] - bbox[0]
    draw.text((center_x - width // 2, y), text, fill=TEXT_COLOR, font=font)


def max_x(node: PathwayNode) -> int:
    """Return the rightmost node x coordinate in the tree."""
    return max([node.x, *(max_x(child) for child in node.children)])


def draw_edges(draw: ImageDraw.ImageDraw, node: PathwayNode) -> None:
    """Draw curved connectors for a pathway tree."""
    for child in node.children:
        start_x, start_y = node.x, node.y
        end_x, end_y = child.x, child.y
        mid_x = (start_x + end_x) // 2
        points = []
        for index in range(31):
            t = index / 30
            u = 1 - t
            curve_x = (
                u**3 * start_x
                + 3 * u**2 * t * mid_x
                + 3 * u * t**2 * mid_x
                + t**3 * end_x
            )
            curve_y = (
                u**3 * start_y
                + 3 * u**2 * t * start_y
                + 3 * u * t**2 * end_y
                + t**3 * end_y
            )
            points.append((int(curve_x), int(curve_y)))

        for first, second in zip(points, points[1:]):
            draw.line([first, second], fill=LINE_COLOR, width=2)
        draw_edges(draw, child)


def draw_node(
    image: Image.Image,
    draw: ImageDraw.ImageDraw,
    node: PathwayNode,
    label_font: ImageFont.ImageFont,
    small_font: ImageFont.ImageFont,
) -> None:
    """Draw one pathway node and its descendants."""
    radius = MAIN_RADIUS if node.is_target else SMALL_RADIUS
    count = max(len(node.molecules), 1)
    top_y = node.y - (count * MOLECULE_CELL_HEIGHT) // 2

    draw_centered_text(draw, node.x, top_y - 24, node.label, label_font)
    for index, molecule_data in enumerate(node.molecules):
        center_y = top_y + index * MOLECULE_CELL_HEIGHT + MOLECULE_CELL_HEIGHT // 2
        fill = TARGET_FILL if node.is_target else STEP_FILL
        stroke = TARGET_STROKE if node.is_target else STEP_STROKE
        draw.ellipse(
            [
                node.x - radius,
                center_y - radius,
                node.x + radius,
                center_y + radius,
            ],
            fill=fill,
            outline=stroke,
            width=2,
        )

        smiles = str(molecule_data.get("smiles", ""))
        molecule_size = int(radius * 1.5)
        molecule_image = render_molecule(smiles, molecule_size + 40)
        if molecule_image is not None:
            molecule_image = molecule_image.resize(
                (molecule_size, molecule_size),
                Image.Resampling.LANCZOS,
            )
            image.paste(
                molecule_image,
                (node.x - molecule_size // 2, center_y - molecule_size // 2),
                molecule_image,
            )

        formula, mass = molecule_metadata(molecule_data)
        draw_centered_text(draw, node.x, center_y + radius + 6, mass, small_font)
        draw_centered_text(draw, node.x, center_y + radius + 22, formula, small_font)

    for child in node.children:
        draw_node(image, draw, child, label_font, small_font)


def visualize_pathway(result: dict[str, Any]) -> Image.Image:
    """Render a formatted pathway dictionary as a PIL image.

    Parameters
    ----------
    result : dict[str, Any]
        Dictionary containing ``steps`` and ``dependencies`` keys.

    Returns
    -------
    PIL.Image.Image
        Rendered pathway image.
    """
    label_font = get_font(15)
    small_font = get_font(13)
    root = build_tree(result)
    if root is None:
        image = Image.new("RGB", (400, 200), BACKGROUND)
        draw = ImageDraw.Draw(image)
        draw_centered_text(draw, 200, 90, "No synthesis steps", label_font)
        return image

    tree_height = node_height(root)
    canvas_height = tree_height + PADDING * 2
    layout(root, PADDING + MAIN_RADIUS, PADDING, PADDING + tree_height)
    canvas_width = max_x(root) + MAIN_RADIUS + PADDING + 60

    image = Image.new("RGBA", (canvas_width, canvas_height), (*BACKGROUND, 255))
    draw = ImageDraw.Draw(image)
    draw_edges(draw, root)
    draw_node(image, draw, root, label_font, small_font)
    return image.convert("RGB")
