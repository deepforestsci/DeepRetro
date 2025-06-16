"""Utility modules for the DeepRetro project.

This package contains various helper modules and functions that support
the core logic of the application. This includes utilities for molecule
manipulation, LLM interactions, stability checks, logging, and more.

Key Modules:
    *   `az`: Utilities related to AiZynthFinder.
    *   `custom_logging`: Configuration for custom logging, including structured
        logging with `structlog` and job-specific file handlers.
    *   `hallucination_checks`: Functions to check for and potentially mitigate
        hallucinations in Large Language Model outputs.
    *   `job_context`: Defines context variables for job-specific information,
        useful for request tracing and logging.
    *   `llm`: Utilities for interacting with various Large Language Models (LLMs),
        such as handling API calls, prompts, and responses.
    *   `parse`: Parsing utilities, for example, for SMILES strings, reaction
        notations, or other chemical data formats.
    *   `stability_checks`: Functions for assessing the chemical stability or
        validity of molecules.
    *   `utils_molecule`: General-purpose utilities for working with molecules,
        likely using RDKit or similar cheminformatics libraries.

To use these utilities, import them as needed, e.g.:
`from src.utils import llm_utils` (if `llm.py` is aliased or contains `llm_utils`)
or `from src.utils.llm import specific_function`.
"""
