#!/usr/bin/env python
"""CLI entry point for the DeepRetro SA/SCScore report generator.

Thin wrapper around :func:`deepretro.report_scores.main`. Example:

    python scripts/generate_score_report.py \\
        --base ~/Downloads/deepretro_batch_<ts> \\
        --molecules molecules.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from deepretro.report_scores import main  # noqa: E402

if __name__ == "__main__":
    main()
