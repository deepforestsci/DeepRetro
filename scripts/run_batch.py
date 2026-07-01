#!/usr/bin/env python
"""CLI entry point for the DeepRetro batch retrosynthesis runner.

Thin wrapper around :func:`deepretro.batch.main` so the pipeline logic stays
importable and unit-tested. Example:

    python scripts/run_batch.py \\
        --sheet-url "https://docs.google.com/spreadsheets/d/<ID>/export?format=csv&gid=<GID>" \\
        --molecules molecules.txt \\
        --out batch_output \\
        --solve-mode single_step_agent \\
        --tool-backend sandbox
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running as a plain script (python scripts/run_batch.py) without an
# installed package: put the repo root on sys.path so ``deepretro`` imports.
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from deepretro.batch import main  # noqa: E402  (after sys.path bootstrap)

if __name__ == "__main__":
    main()
