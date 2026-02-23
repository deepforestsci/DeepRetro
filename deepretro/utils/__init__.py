"""Utilities for retrosynthesis: LLM calls, AiZynthFinder, and feature extraction.

Submodules
----------
az : AiZynthFinder integration (run_az, run_az_with_img, is_basic_molecule)
feat_extract : Domain feature extraction for reaction-step featurization
llm : LLM calls, response parsing, and retrosynthesis pipeline
"""

from deepretro.utils.az import is_basic_molecule, run_az, run_az_with_img
from deepretro.utils.feat_extract import (
    NUM_DOMAIN_FEATURES,
    extract_domain_features_single,
)
from deepretro.utils.llm import (
    build_addon_prompt,
    build_completion_params,
    call_LLM,
    extract_tag_content,
    llm_pipeline,
    obtain_prompt,
    parse_response,
    validate_json_response,
)

__all__ = [
    "build_addon_prompt",
    "build_completion_params",
    "call_LLM",
    "extract_domain_features_single",
    "extract_tag_content",
    "is_basic_molecule",
    "llm_pipeline",
    "NUM_DOMAIN_FEATURES",
    "obtain_prompt",
    "parse_response",
    "run_az",
    "run_az_with_img",
    "validate_json_response",
]
