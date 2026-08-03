# Dev-Added Package Tests and Docstrings Design

## Objective

Improve test coverage and API documentation for production Python files newly
added on the DeepRetro `dev` branch. The cleanup will add meaningful tests and
NumPy-style docstrings with realistic examples. Production-file edits are
limited strictly to docstring text and must not change executable code.

## Frozen Diff Boundary

The eligibility boundary is the three-dot diff between these refreshed remote
references on 2026-08-03:

- `upstream/main`: `dc8fcad`
- `upstream/dev`: `0575a1f`
- merge base: `dc8fcad`

Eligible production files are the `.py` files reported by:

```powershell
git diff --name-only --diff-filter=A upstream/main...upstream/dev -- deepretro
```

This yields:

1. `deepretro/agents/__init__.py`
2. `deepretro/agents/loop.py`
3. `deepretro/agents/message_history.py`
4. `deepretro/agents/sandbox.py`
5. `deepretro/agents/tools.py`
6. `deepretro/algorithms/autosolve.py`
7. `deepretro/batch.py`
8. `deepretro/models/hallucination_checker.py`
9. `deepretro/models/hallucination_trainer.py`
10. `deepretro/models/hallucination_utils.py`
11. `deepretro/report_scores.py`
12. `deepretro/score.py`

Files that existed on `main` and were only modified on `dev` are excluded.
Examples include `deepretro/metadata.py`, `deepretro/utils/llm_helpers.py`, and
`deepretro/featurizers/reactionstep.py`.

Tests and project documentation may be added or edited as needed to verify and
document the eligible production files. Unrelated dirty-worktree files are out
of scope.

## Allowed Change Types

The final implementation diff may contain only:

- docstring additions or revisions inside eligible production files;
- test additions or revisions under `deepretro/tests/`; and
- project documentation for this cleanup.

Executable statements, imports, annotations, decorators, signatures, constants,
formatting outside docstrings, configuration, and runtime behavior must remain
unchanged. If a new test exposes a production defect, the defect will be
reported separately instead of fixed in this cleanup.

## Coverage Boundary

The cleanup covers:

- every public class;
- every public method;
- every public function; and
- private helpers only when they contain meaningful branching, parsing,
  aggregation, resource handling, or error handling.

Tests must assert externally meaningful output, state, or failure behavior.
Tests that merely invoke a function for line coverage are not acceptable.
Already-strong agent and autosolve coverage will not be duplicated.

## Documentation Contract

Eligible public APIs will have clear NumPy-style docstrings that describe their
purpose, parameters, returns, raised exceptions where relevant, and behavioral
constraints. Each public API will include a realistic `Examples` section.

Examples should be executable wherever practical. They must use deterministic
local inputs and must not require provider credentials, network access, live
LLM calls, Langfuse, or a live AiZynthFinder installation. Where executing the
full integration is inherently external or expensive, the example will show a
valid local construction or data-transformation path without pretending to run
the unavailable integration.

Docstrings must describe existing behavior. This task will not alter production
semantics or non-docstring code to make an example easier to write.

## Test Strategy

### Scoring

`deepretro/score.py` is the highest-priority gap because it has no committed
package test module on `dev`. Tests will cover the public scoring API,
canonicalization, malformed inputs, optional-scorer fallbacks, step/pathway
summaries, and stable empty-result schemas using real local chemistry behavior
where available.

### Hallucination model modules

`deepretro/models/hallucination_checker.py`,
`deepretro/models/hallucination_trainer.py`, and
`deepretro/models/hallucination_utils.py` will receive direct public-contract
coverage. Tests may use small real local datasets, temporary files, and locally
available model objects. Expensive training runs and external integrations are
outside the unit-test boundary.

### Batch and score reports

`deepretro/batch.py` and `deepretro/report_scores.py` will use real temporary
CSV, JSON, Markdown, and LaTeX files to cover parsing, aggregation, rendering,
output naming, and deterministic failure behavior.

### Agents and autosolve

The agent, sandbox, tool-registry, message-history, and autosolve modules already
have substantial direct coverage. They will be audited for public-contract gaps,
but new tests will be added only for meaningful missing behavior. Their examples
will be normalized where needed without duplicating existing tests.

### Lazy exports

`deepretro/agents/__init__.py` will be documented and tested only to the extent
that its lazy exports are part of the supported public package contract.

## No-Mock and Integration Rules

New tests must not use:

- `unittest.mock`, `Mock`, or `MagicMock`;
- `patch` or `pytest.monkeypatch`;
- simulated provider clients; or
- live LLM, cloud, Langfuse, or AiZynthFinder calls.

Tests will instead exercise pure logic, validation boundaries, real local
objects, and temporary filesystem inputs. External behavior that cannot be
tested under those constraints will remain outside this unit-test cleanup.

## Workflow

1. Create a clean isolated worktree from the refreshed `upstream/dev` state.
2. Run the existing package suite to establish a baseline.
3. Stop and report exact failures if the clean baseline does not pass.
4. Audit eligible public APIs and behavior-heavy private helpers.
5. Add each missing behavior test before changing its associated documentation.
6. Add or refine NumPy-style docstrings and examples.
7. Run focused validation after each module group.
8. Run the full package validation and final scope audit.

No commit will contain unrelated working-tree changes. Nothing will be pushed
and no pull request will be created without Kamal's confirmation.

## Verification

Verification will include:

- focused tests for each changed module group;
- the complete `deepretro/tests` suite with a repo-local `--basetemp`;
- practical doctest/example execution;
- the package's configured Ruff checks;
- the package's configured type checks;
- `python -m py_compile` for changed Python files;
- `git diff --check`;
- a scope audit showing production edits are limited to docstrings in eligible
  files;
- a scan proving newly added tests contain no forbidden mocking constructs.

## Success Criteria

The cleanup is complete when:

1. all eligible public APIs have accurate NumPy-style docstrings and examples;
2. meaningful public behavior and important private branching are directly
   tested without mocks or live services;
3. production-file changes contain docstrings only, with no executable-code,
   configuration, formatting, or behavior changes;
4. focused and package-wide validation pass, or any environment limitation is
   isolated and reported with exact evidence;
5. the final diff contains only approved production files, corresponding tests,
   and project documentation; and
6. all existing unrelated user changes remain untouched.
