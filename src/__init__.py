"""Main source package for the DeepRetro application.

This package contains the core logic, API endpoints, processing pipelines,
and utility modules that constitute the DeepRetro application.

Key submodules and packages include:
    api: Flask application providing the HTTP API for retrosynthesis.
    main: Core orchestration logic for running retrosynthesis tasks.
    utils: A sub-package containing various helper utilities for cheminformatics,
           caching, logging, etc.
    runners: Scripts for batch processing and running experiments (not a package itself).
    variables: Definitions of global constants, mappings, and predefined data.
    cache: Caching functionalities.
    metadata: Functions related to metadata generation or processing.
    prithvi / rec_prithvi: Modules related to specific models or processing flows,
                         likely involving LLM interactions.

This `__init__.py` makes `src` a Python package.
"""
