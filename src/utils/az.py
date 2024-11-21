"""Utils for AiZynthFinder"""
import os
from aizynthfinder.aizynthfinder import AiZynthFinder
from typing import Any, Dict, List, Optional, Sequence
from src.variables import BASIC_MOLECULES, ENCODING_SCALABILITY
from src.cache import cache_results
import rootutils

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

AZ_MODEL_CONFIG_PATH = f"{root_dir}/{os.getenv('AZ_MODEL_CONFIG_PATH')}"
finder = AiZynthFinder(configfile=AZ_MODEL_CONFIG_PATH)
finder.stock.select("zinc")
finder.expansion_policy.select("uspto")
finder.filter_policy.select("uspto")


@cache_results
def run_az(smiles: str) -> tuple[Any, Sequence[Dict[str, Any]]]:
    """Run the retrosynthesis using AiZynthFinder

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule

    Returns
    -------
    tuple[Any, Sequence[Dict[str, Any]]]
        A tuple containing the status of the retrosynthesis, 
        the results dictionary
    """
    # if simple molecule, skip the retrosynthesis
    if smiles in BASIC_MOLECULES:
        return True, [{
            'type': 'mol',
            'hide': False,
            'smiles': smiles,
            'is_chemical': True,
            'in_stock': True,
        }]
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats = finder.extract_statistics()
    status = stats['is_solved']
    result_dict = finder.routes.dict_with_extra(include_metadata=True,
                                                include_scores=True)
    return status, result_dict
