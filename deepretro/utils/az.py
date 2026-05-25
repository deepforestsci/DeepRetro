"""AiZynthFinder integration for template-based retrosynthesis.

Runs AiZynthFinder on target molecules, with optional image export.
Uses ZINC stock and USPTO expansion/filter policies by default.
Requires ``AZ_MODEL_CONFIG_PATH`` or ``AZ_MODELS_PATH`` environment variables.
"""

import os
from pathlib import Path
from typing import Any, Dict, Sequence

from dotenv import load_dotenv
from PIL.Image import Image
from rdkit import Chem
from rdkit.Chem import rdqueries

from deepretro.utils.variables import BASIC_MOLECULES

load_dotenv()


def _find_project_root(start: Path, marker: str = ".project-root") -> Path:
    """Walk up from *start* until a directory containing *marker* is found."""
    for directory in [start.resolve(), *start.resolve().parents]:
        if (directory / marker).exists():
            return directory
    return Path.cwd()


root_dir = _find_project_root(Path(__file__))

ENABLE_LOGGING = (
    False if os.getenv("ENABLE_LOGGING", "true").lower() == "false" else True
)

# Paths from env; required for AiZynthFinder config and model files
AZ_MODEL_CONFIG_PATH = f"{root_dir}/{os.getenv('AZ_MODEL_CONFIG_PATH')}"
AZ_MODELS_PATH = f"{root_dir}/{os.getenv('AZ_MODELS_PATH')}"


def _log(message: str, logger=None):
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


def run_az(
    smiles: str, az_model: str = "USPTO"
) -> tuple[bool, Sequence[Dict[str, Any]]]:
    """Run the retrosynthesis using AiZynthFinder.

    Example
    -------
    >>> from deepretro.utils.az import run_az
    >>> status, result_dict = run_az("C1CCCCC1", "USPTO")  # doctest: +SKIP
    >>> isinstance(status, bool) and isinstance(result_dict, list)  # doctest: +SKIP
    True

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.
    az_model : str, optional
        AiZynthFinder model variant (e.g. ``"USPTO"``, ``"Pistachio_50"``),
        by default ``"USPTO"``.

    Returns
    -------
    tuple[bool, Sequence[Dict[str, Any]]]
        ``(solved, routes)`` — whether a route was found and the route data.
    """
    try:
        config_path = f"{AZ_MODELS_PATH}/{az_model}/config.yml"
        with open(config_path, "r") as _:
            config_filename = config_path
    except FileNotFoundError:
        _log(f"AZ_MODEL_CONFIG_PATH not found at {config_path}")
        try:
            with open(AZ_MODEL_CONFIG_PATH, "r") as _:
                config_filename = AZ_MODEL_CONFIG_PATH
        except FileNotFoundError:
            raise FileNotFoundError(
                f"AZ_MODEL_CONFIG_PATH not found at {AZ_MODEL_CONFIG_PATH}"
            )
    # if simple molecule, skip the retrosynthesis
    if smiles in BASIC_MOLECULES or is_basic_molecule(smiles):
        return True, [
            {
                "type": "mol",
                "hide": False,
                "smiles": smiles,
                "is_chemical": True,
                "in_stock": True,
            }
        ]
    from aizynthfinder.aizynthfinder import AiZynthFinder

    finder = AiZynthFinder(configfile=config_filename)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats = finder.extract_statistics()
    status = bool(stats["is_solved"])
    result_dict = finder.routes.dict_with_extra(
        include_metadata=True, include_scores=True
    )
    return status, result_dict


def run_az_with_img(
    smiles: str,
) -> tuple[bool, Sequence[Dict[str, Any]], Sequence[Image | None] | None]:
    """Run the retrosynthesis using AiZynthFinder.

    Example
    -------
    >>> from deepretro.utils.az import run_az_with_img
    >>> status, result_dict, images = run_az_with_img("C1CCCCC1")  # doctest: +SKIP
    >>> isinstance(status, bool)  # doctest: +SKIP
    True

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.

    Returns
    -------
    tuple[bool, Sequence[Dict[str, Any]], Sequence[Image] | None]
        ``(solved, routes, images)`` — solved status, route data, and
        optional route images (PNG bytes). Uses ``AZ_MODEL_CONFIG_PATH``.
    """
    # if simple molecule, skip the retrosynthesis
    if smiles in BASIC_MOLECULES or is_basic_molecule(smiles):
        return (
            True,
            [
                {
                    "type": "mol",
                    "hide": False,
                    "smiles": smiles,
                    "is_chemical": True,
                    "in_stock": True,
                }
            ],
            None,
        )
    from aizynthfinder.aizynthfinder import AiZynthFinder

    finder = AiZynthFinder(configfile=AZ_MODEL_CONFIG_PATH)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats: dict[str, Any] = finder.extract_statistics()
    status = bool(stats["is_solved"])
    result_dict: Sequence[dict[str, Any]] = finder.routes.dict_with_extra(
        include_metadata=True, include_scores=True
    )
    images: Sequence[Image | None] = finder.routes.images
    return status, result_dict, images


def is_basic_molecule(smiles: str) -> bool:
    """Check if the molecule is a basic molecule
    (if number of C atoms is less than 5 or total atoms < 5).

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule

    Returns
    -------
    bool
        True if the molecule is a basic molecule, False otherwise

    Examples
    --------
    >>> from deepretro.utils.az import is_basic_molecule
    >>> is_basic_molecule("C")
    True
    >>> is_basic_molecule("CC")
    True
    >>> is_basic_molecule("C1CCCCC1")
    False
    >>> is_basic_molecule("CCO")
    True
    >>> is_basic_molecule("invalid_smiles!!")
    False
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return False
    if mol is None:
        return False

    q = rdqueries.AtomNumEqualsQueryAtom(6)
    num_c_atoms = len(mol.GetAtomsMatchingQuery(q))
    # if total number of atoms is less than 5, return True
    if mol.GetNumAtoms() < 5:
        return True
    elif num_c_atoms < 5:
        return True
    # if total number of C atoms is less than 5, return True

    return False
