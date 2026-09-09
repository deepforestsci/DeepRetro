# DeepRetro

DeepRetro provides retrosynthesis search, reaction featurizers, and heuristic
reaction checks. The heuristic checker reports structural warnings and ranks
candidate reactions; its score is not a calibrated probability of chemical
correctness.

## Installation

Add the package to a Python project using uv:

```sh
uv add /path/to/DeepRetro/deepretro
```

## Heuristic scoring

This local example requires no API credentials or trained model:

```python
from deepretro.algorithms.hallucination_checker import calculate_hallucination_score

report = calculate_hallucination_score(reactant_smiles="CCO", product_smiles="CC=O")
print(report["score"], report["severity"])
```

Search can retain flagged candidates as fallbacks, so inspect reaction warnings
before interpreting a completed route. Heuristic checks cannot establish
experimental feasibility or assess unspecified reaction conditions.

See the [project documentation](https://github.com/deepforestsci/DeepRetro/tree/dev/docs)
for solver configuration, custom weights, and development instructions.
