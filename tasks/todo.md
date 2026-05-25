# AutoSolver Package Port

## Plan

- [x] Create a new branch from `porting/metadata-llm`.
- [x] Inspect `origin/autosolve` and identify stale or duplicated migration code.
- [x] Port the AutoSolver public API into the `deepretro` package layout.
- [x] Replace `src.*` imports and duplicated utilities with package-local modules.
- [x] Clean implementation style to match current typed, lazy-import package standards.
- [x] Add focused unit tests for solver behavior, validation, recursion, and image wiring.
- [x] Run formatter, linter, type checker where practical, and targeted tests.
- [x] Ask Claude Code to review the tests with `claude -p`.
- [x] Push the branch and open a draft PR.

## Success Criteria

- `deepretro.AutoSolver` and `deepretro.algorithms.AutoSolver` import lazily.
- AutoSolver can use AiZynthFinder results directly and LLM fallback when AZ fails.
- Cycles, invalid SMILES, max-depth exits, and all-unsolved pathways are deterministic.
- Tests do not require network calls, real LLM credentials, or AiZynthFinder model files.
- The PR is draft and describes validation evidence plus known gaps.

## Review

- Ported AutoSolver without protecting-group prompt support; that scope is
  reserved for a later PR.
- Used the existing `deepretro.utils.llm.llm_pipeline` instead of adding the
  stale duplicate `deepretro.algorithms.llm` module from `origin/autosolve`.
- Added focused AutoSolver and hallucination helper tests, then expanded them
  after Claude review to cover partial fallback trees and multi-reactant
  pathways.
- Validation passed:
  - `PYTHONPATH=. uv run --project deepretro --extra dev python -m pytest deepretro/tests/test_autosolve.py deepretro/tests/test_hallucination_helpers.py`
  - `PYTHONPATH=. uv run --project deepretro --extra dev python -m pytest deepretro/tests -m "not slow"`
  - `PYTHONPATH=. uv run --project deepretro --extra dev ruff format --check ...`
  - `PYTHONPATH=. uv run --project deepretro --extra dev ruff check ...`
  - `PYTHONPATH=. uv run --project deepretro --extra dev ty check deepretro/algorithms/autosolve.py deepretro/models/hallucination_helpers.py deepretro/utils/visualize.py`
  - `PYTHONPATH=. uv run --project deepretro --extra dev --extra docs sphinx-build -b html docs/source docs/_build/html`
