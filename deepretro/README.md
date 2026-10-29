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

## Autosolve metadata

Parsed routes report completion separately from reaction provenance:

```python
result = solver.autosolve(target_smiles)  # configured AutoSolver instance
print(result["solved"])  # AZ succeeded for the target or every terminal branch
print(result["az_solved"])  # complete route with only AZ reactions
print([step["solved_by"] for step in result["steps"]])  # "llm" or "az"
```

A completed mixed LLM/AZ route has `solved=True` and `az_solved=False`.
`az_summary.az_solved_all` means AZ closed all terminal branches, even if
some reactions came from the LLM. AZ stock hits have no reaction steps.

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

## Guarded heuristic search

`GuardedSolver` checks both template-generated and LLM reactions, rejects flagged
fallbacks, and sends rejection reasons into bounded repair calls. It requires an
explicit stock lookup and writes partial-route checkpoints and repair events.
This offline example uses an illustrative stock and proposal runner:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from deepretro import AutoSolver
from deepretro.algorithms.guarded import GuardedSolver

backend = AutoSolver(
    hallucination_mode="heuristic",
    az_runner=lambda smiles, model: (False, []),
    llm_runner=lambda molecule, **kwargs: ([["CCO"]], ["oxidation"], [0.9]),
)
with TemporaryDirectory() as directory:
    solver = GuardedSolver(backend, {"CCO"}.__contains__, Path(directory) / "run")
    tree, passed = solver.solve("CC=O")
    result = backend.parse(tree, solved=passed)
```

For production, supply a catalog-backed callback instead of the example set.
`deepretro.utils.disk_stock.SortedInchiKeyQuery` provides exact, memory-mapped
InChIKey lookup for a producer-verified sorted `.npy` stock. Heuristic acceptance
does not establish chemical feasibility. Ordinary `AutoSolver` behavior is
unchanged; guarded search currently supports pipeline mode only.

### Agent feedback

The agent's `check_hallucination` tool returns heuristic `score`, `severity`,
`penalties`, and `explanation.detected_issues`. Inspect the same assessment locally:

```python
from deepretro.models.hallucination_checker import HallucinationChecker

checker = HallucinationChecker(checker_type="heuristic")
report = checker.assess("c1ccccc1", ["CC"])
print(report["flagged"])  # True, even if search retains this fallback
print(report["explanation"]["detected_issues"])
```

The tool exposes `flagged` as `is_hallucination`. For unassessable inputs it is
`None` (`null` in JSON); inspect `message` for the reason. Explanations include
ring-size and relative substituent-position changes, without exact atom/bond
localization.


## Run package tests

From the repository root:

```sh
uv sync --project deepretro --extra dev --locked
uv run --project deepretro --no-sync pytest deepretro/tests
```

Unit tests block Python socket connections in the test process. Integration
tests are skipped even when provider credentials are present. External programs
and child processes still require explicit test doubles.

To request live provider or model tests deliberately:

```sh
uv run --project deepretro --no-sync pytest deepretro/tests --run-integration -m integration
```

Provider tests require the appropriate credentials and may incur charges.
AiZynthFinder tests require the `az` extra and an explicitly provisioned model
directory in `AZ_TEST_MODELS_DIR`, containing `config.yml`. Test fixtures never
download models. A missing setting skips these tests; a configured directory
without `config.yml` fails with an actionable error.
