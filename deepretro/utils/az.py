"""AiZynthFinder integration for template-based retrosynthesis.

Runs AiZynthFinder on target molecules, with optional image export.
Uses ZINC stock and USPTO expansion/filter policies by default.
Requires ``AZ_MODEL_CONFIG_PATH`` or ``AZ_MODELS_PATH`` environment variables.
Caching is opt-in through an explicit ``CacheManager`` argument.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Sequence, cast

import structlog
from rdkit import Chem
from rdkit.Chem import rdqueries

from deepretro.utils.cache import CacheManager, make_cache_key
from deepretro.utils.variables import BASIC_MOLECULES

if TYPE_CHECKING:
    from PIL.Image import Image

try:
    from aizynthfinder.aizynthfinder import AiZynthFinder
except ImportError:
    AiZynthFinder = None  # type: ignore[assignment, misc]

logger = structlog.get_logger()

PROJECT_ROOT = Path(__file__).resolve().parents[2]

AZ_MODEL_CONFIG_PATH = f"{PROJECT_ROOT}/{os.getenv('AZ_MODEL_CONFIG_PATH')}"
AZ_MODELS_PATH = f"{PROJECT_ROOT}/{os.getenv('AZ_MODELS_PATH')}"


def _basic_molecule_route(smiles: str) -> list[Dict[str, Any]]:
    """Return the solved route payload for a feedstock/basic molecule."""
    return [
        {
            "type": "mol",
            "hide": False,
            "smiles": smiles,
            "is_chemical": True,
            "in_stock": True,
        }
    ]


def _resolve_config(az_model: str | None = None) -> str:
    """Resolve AiZynthFinder config path, falling back to ``AZ_MODEL_CONFIG_PATH``."""
    if az_model is not None:
        config_path = f"{AZ_MODELS_PATH}/{az_model}/config.yml"
        try:
            with open(config_path, "r") as _:
                return config_path
        except FileNotFoundError:
            logger.warning("AZ config not found, trying fallback", path=config_path)

    try:
        with open(AZ_MODEL_CONFIG_PATH, "r") as _:
            return AZ_MODEL_CONFIG_PATH
    except FileNotFoundError:
        raise FileNotFoundError(
            f"AZ_MODEL_CONFIG_PATH not found at {AZ_MODEL_CONFIG_PATH}"
        )


def _run_finder(config_filename: str, smiles: str) -> Any:
    """Set up and run AiZynthFinder, returning the finder instance."""
    if AiZynthFinder is None:
        raise ImportError(
            "AiZynthFinder support requires optional dependencies. "
            "Install the package with `deepretro[az]`."
        )
    finder = AiZynthFinder(configfile=config_filename)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    return finder


def run_az(
    smiles: str,
    az_model: str = "USPTO",
    cache: CacheManager | None = None,
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
    cache : CacheManager | None, optional
        Explicit cache instance used to memoize results for this call. When
        ``None``, no cache is read or written.

    Returns
    -------
    tuple[bool, Sequence[Dict[str, Any]]]
        ``(solved, routes)`` — whether a route was found and the route data.

    Notes
    -----
    Install the package with ``deepretro[az]``. Caching is disabled unless an
    explicit ``cache=CacheManager(...)`` is supplied.
    """
    cache_key = make_cache_key("run_az", smiles, az_model=az_model, version=1)
    cache_miss = object()
    if cache is not None:
        cached_result = cache.get(cache_key, default=cache_miss)
        if cached_result is not cache_miss:
            return cast(tuple[bool, Sequence[Dict[str, Any]]], cached_result)

    if smiles in BASIC_MOLECULES or is_basic_molecule(smiles):
        result = (True, _basic_molecule_route(smiles))
        if cache is not None:
            cache.set(cache_key, result, tag=smiles)
        return result

    config_filename = _resolve_config(az_model)
    finder = _run_finder(config_filename, smiles)
    stats = finder.extract_statistics()
    status = bool(stats["is_solved"])
    result_dict = finder.routes.dict_with_extra(
        include_metadata=True, include_scores=True
    )
    result = (status, result_dict)
    if cache is not None:
        cache.set(cache_key, result, tag=smiles)
    return result


def run_az_with_img(
    smiles: str,
    cache: CacheManager | None = None,
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
    cache : CacheManager | None, optional
        Explicit cache instance used to memoize results for this call. When
        ``None``, no cache is read or written.

    Returns
    -------
    tuple[bool, Sequence[Dict[str, Any]], Sequence[Image] | None]
        ``(solved, routes, images)`` — solved status, route data, and
        optional route images (PNG bytes). Uses ``AZ_MODEL_CONFIG_PATH``.

    Notes
    -----
    Install the package with ``deepretro[az]``. Caching is disabled unless an
    explicit ``cache=CacheManager(...)`` is supplied.
    """
    cache_key = make_cache_key("run_az_with_img", smiles, version=1)
    cache_miss = object()
    if cache is not None:
        cached_result = cache.get(cache_key, default=cache_miss)
        if cached_result is not cache_miss:
            return cast(
                tuple[bool, Sequence[Dict[str, Any]], Any],
                cached_result,
            )

    if smiles in BASIC_MOLECULES or is_basic_molecule(smiles):
        result = (True, _basic_molecule_route(smiles), None)
        if cache is not None:
            cache.set(cache_key, result, tag=smiles)
        return result

    config_filename = _resolve_config()
    finder = _run_finder(config_filename, smiles)
    stats: dict[str, Any] = finder.extract_statistics()
    status = bool(stats["is_solved"])
    result_dict: Sequence[dict[str, Any]] = finder.routes.dict_with_extra(
        include_metadata=True, include_scores=True
    )
    images: Sequence[Image | None] = finder.routes.images
    result = (status, result_dict, images)
    if cache is not None:
        cache.set(cache_key, result, tag=smiles)
    return result


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
