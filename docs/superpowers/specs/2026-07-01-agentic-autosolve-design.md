# Agentic tool-calling retrosynthesis in `deepretro/`

**Status:** Approved design — pending 4-persona review before implementation
**Branch:** `feat/agentic-autosolve` (off `dev`) → PR to `dev`
**Date:** 2026-07-01

## 1. Goal

Adapt DeepRetro's recursive retrosynthesis algorithm to run as an **agentic loop with
tool calling**, living inside the `deepretro/` package. Three deliverables:

1. Port the recursive loop + metadata into `deepretro/` (built on the unmerged
   `origin/autosolve-pipeline` branch, rewired to current `dev` APIs).
2. An agentic layer where the LLM can call tools — **molecule validity**, **stability**,
   **hallucination** (template), and **write-and-run-Python-in-a-sandbox** — with tool
   calling and sandbox code-execution both selectable.
3. A batch script: download a CSV from Google Sheets → train the hallucination checker
   (template) → run molecules from a `.txt` file through the full algorithm → dump
   outputs as `folder/<timestamp>/<molecule>/pathway_<i>.json`.

## 2. Current state (verified against `dev`)

Already in `deepretro/` on `dev`:

- **LLM layer** (100% `litellm`): `utils/llm_interface.py` (provider classes + `call()`),
  `utils/llm.py` (`call_LLM`, `llm_pipeline` — single-step retrosynthesis with rising-temperature
  retries + validity/stability/hallucination filters), `utils/llm_helpers.py`
  (`build_completion_params` — the one place that assembles `litellm.completion(**params)`).
- **Checks**: `algorithms/stability_checker.py::check_molecule_stability`,
  `algorithms/hallucination_checker.py::calculate_hallucination_score` (heuristic),
  `utils/utils_molecule.py::is_valid_smiles` / `validity_check`,
  `algorithms/pipeline_checks.py::hallucination_checker` (pipeline hook).
- **ML hallucination**: `models/hallucination_classifier.py::HallucinationClassifier`
  (XGBoost via DeepChem `GBDTModel`; `fit/evaluate/predict_probability/save/load`),
  `models/hallucination_helpers.py::MLChecker` + `resolve_hallucination(mode, classifier)`.
- **Metadata**: `metadata.py::recommend_reaction_metadata` (reagent→conditions→literature agents).
- **AZ**: `utils/az.py::run_az(smiles, az_model)`. **Parse**: `utils/parse.py::format_output`.
- **Training data path**: `featurizers/reactionstep.py::ReactionStepFeaturizer`,
  `data/loader.py::ReactionDataLoader` (cols `product`/`reactants`/`label`) + `stratified_split`.

Missing on `dev`: the recursive loop, any tool-calling/agent layer, any sandbox. There is a
dead `TOOLS` list in `utils/variables.py` (never imported) and a one-tool prototype on
`origin/agentic-solution` (`agentic_call_LLM` + `src/utils/tools.py`).

## 3. Reuse decision: build on `origin/autosolve-pipeline`

That branch's `deepretro/algorithms/autosolve.py` (~636 lines) already implements the
`AutoSolver` (recursive `solve`, `single_step`, `parse`, `add_metadata`, `autosolve`,
`run_llm`) with dependency injection and ~660 lines of tests. It is **not** cleanly rebasable —
it targets the older `autosolve-foundations` API that `dev`'s merged hallucination layer
replaced. Reuse it as the structural basis with a **bounded adaptation layer**:

| Branch import | Status on `dev` | Adaptation |
|---|---|---|
| `hallucination_helpers.filter_with_checker` | absent (replaced by `MLChecker`) | call `checker(product, pathways) -> (status, kept)` where `checker = resolve_hallucination(mode, clf)` |
| `hallucination_helpers.HallucinationPredictor` | absent | type as the checker `Callable`; drop the alias |
| `utils.typing.HallucinationChecker` | absent | add alias `Callable[[str, list], tuple[int, list]]` |
| `utils_molecule.canonicalize` | absent | port the trivial RDKit wrapper from foundations |
| `run_az`, `llm_pipeline`, `format_output`, `recommend_reaction_metadata`, recommender types | present | no change |

## 4. Architecture

```
deepretro/
├── algorithms/autosolve.py     port AutoSolver (rewired) + reasoning_mode / tool_backend flags
├── agents/                     NEW subpackage
│   ├── __init__.py
│   ├── tools.py                tool registry: schemas + executors + dispatcher
│   ├── sandbox.py              subprocess sandbox behind a Sandbox protocol
│   └── loop.py                 agentic_single_step() (full) + agentic_orchestrator() (scaffold)
├── utils/utils_molecule.py     + canonicalize()
└── utils/typing.py             + HallucinationChecker alias
scripts/run_batch.py            NEW batch driver
```

### 4.1 `AutoSolver` flags (item: "both, flag-selectable")

`AutoSolver.__init__` gains:

- `reasoning_mode: {"pipeline", "single_step_agent", "orchestrator"} = "pipeline"`
- `tool_backend: {"tools", "sandbox"} = "tools"`
- `sandbox: Sandbox | None = None` (injectable; defaults to the subprocess backend)

Behavior:

- **`pipeline`** (default) — unchanged: `run_llm` calls `llm_pipeline`. Nothing is ripped out.
- **`single_step_agent`** (fully built) — `run_llm`/`single_step` calls
  `deepretro.agents.loop.agentic_single_step`, which lets the LLM call the tools to
  self-correct its proposed precursors, then returns `(pathways, explanations, confidence)`
  in the **same shape** `llm_pipeline` returns. The deterministic validity/stability/
  hallucination filters still run afterward as a safety net — the agent's self-checks
  augment, never replace, them.
- **`orchestrator`** (scaffold, flag-wired, minimal impl) — `deepretro.agents.loop.agentic_orchestrator`
  gives one agent tool-driven control of the whole search. Interface + wiring present;
  `single_step_agent` is the deep implementation.

The recursion in `solve()` is unchanged regardless of mode (§5).

### 4.2 Tools (`agents/tools.py`)

Ported/extended from `origin/agentic-solution:src/utils/tools.py`: Anthropic/litellm tool
schema dicts + executor functions + a `TOOL_DEFINITIONS`/`TOOL_EXECUTORS` registry +
`execute_tool(name, input) -> dict` dispatcher + `get_available_tools()`.

| tool | signature | wraps | returns |
|---|---|---|---|
| `validate_smiles` | `(smiles, context?)` | `is_valid_smiles` + RDKit canonical | `{valid, canonical_smiles, error, context}` |
| `check_stability` | `(smiles)` | `check_molecule_stability` | `{assessment, stability_score, issues}` |
| `check_hallucination` | `(product, reactants)` | `resolve_hallucination`/`MLChecker` | **TEMPLATE**: `{configured: false, fallback_severity, note}` when no classifier; real verdict when one is injected. Flagged for the other contributor. |
| `run_python` | `(code)` | `agents/sandbox.py` | `{ok, stdout, stderr, error}` — RDKit + `deepretro.algorithms` checkers importable |

The tool registry is built once and passed into the agent loop; `check_hallucination`'s
executor is bound to whatever `resolve_hallucination(mode, classifier)` returns so the ML
classifier can be swapped in later without touching tool code.

### 4.3 Sandbox (`agents/sandbox.py`) — subprocess backend

```python
class Sandbox(Protocol):
    def run(self, code: str, *, timeout_s: float, mem_mb: int) -> SandboxResult: ...

class SubprocessSandbox:  # default
    # writes code to a temp file, runs `python <file>` in a subprocess with:
    #   - resource.setrlimit(RLIMIT_AS, mem_mb) and RLIMIT_CPU (POSIX)
    #   - a fresh temp cwd, cleaned up after
    #   - captured stdout/stderr, hard wall-clock timeout (kill on expiry)
    # returns SandboxResult(ok, stdout, stderr, error)
```

`ponytail:` comment marks the isolation ceiling — rlimit + no-inherited-fds is adequate for
**LLM-generated, non-adversarial** chem code, but is **not** namespace/network isolation;
upgrade path is a `PodmanSandbox` implementing the same `Sandbox` protocol if code ever
becomes untrusted. Network is best-effort suppressed (cleared proxy env; documented ceiling).

### 4.4 Agent loop (`agents/loop.py`)

Ported/extended from `origin/agentic-solution:agentic_call_LLM`. Uses litellm's OpenAI-style
tool calling through the existing `build_completion_params` seam:

```python
def agentic_single_step(
    molecule, model, *, tools, tool_backend, max_iterations=6,
    llm_runner=None, enable_thinking=True, max_output_tokens=None,
) -> tuple[list[Pathway], list[str], list[float]]:
    # loop: completion(**params, tools=tool_defs)
    #   while finish_reason == "tool_use" / message.tool_calls and iter < max_iterations:
    #       for each tool_call: execute_tool(...) -> append tool result message
    #   on final answer: parse into (pathways, explanations, confidence)  [same shape as llm_pipeline]
```

`tool_backend="sandbox"` restricts the exposed tools to `run_python` (+ read-only checks);
`tool_backend="tools"` exposes the structured check tools. `llm_runner` is injectable for tests
(a fake returning a scripted `tool_use` then a final answer — no network). Iteration cap and
per-tool error capture prevent runaway loops.

## 5. Data flow (recursion terminates when AZ solves)

```
autosolve(smiles):
    route_tree, solved = solve(smiles, depth=0, visited=set())
    output = parse(route_tree)                 # format_output
    output = add_metadata(output)              # recommend_reaction_metadata
    return output, solved

solve(molecule, depth, visited):
    m = canonicalize(molecule)
    if depth >= max_depth or m in visited:  return unsolved_leaf(m), False   # safety guards only
    visited.add(m)
    solved, routes = run_az(m, az_model)
    if solved:  return az_route(routes), True          # ✅ AZ found building blocks → branch ends
    # AZ could not solve → go deeper
    pathways, expl, conf = single_step(m)              # pipeline OR agent, per reasoning_mode
    pathways = validity / stability / hallucination filters
    for pathway in pathways (best-first):
        child_results = [solve(p, depth+1, visited) for p in pathway]
        if all(child solved):  return assemble(pathway, child_results), True
    return best_effort_partial, False
```

A branch stops **only** when AZ solves the node (purchasable/building-block leaf). A full
pathway is solved iff every leaf bottoms out at an AZ-solved node. `max_depth`/cycle-detection
are guards against infinite recursion when AZ never solves a path. This is exactly
`src/rec_prithvi.py::rec_run_prithvi`'s behavior, preserved by `autosolve-pipeline::solve`.

## 6. Batch script (`scripts/run_batch.py`)

CLI:
```
python scripts/run_batch.py \
  --sheet-url <google-sheets /export?format=csv URL>   # REQUIRED (public export, no auth)
  --molecules <path/to/molecules.txt>                  # REQUIRED (one SMILES per line)
  --out <output-dir> \
  --reasoning-mode {pipeline,single_step_agent,orchestrator} \
  --tool-backend {tools,sandbox} \
  --model <litellm model id> \
  --classifier <optional path to trained HallucinationClassifier> \
  --top-k <int, default 3>
```

Steps:

1. `download_sheet_csv(export_url, dest)` — plain `httpx`/`requests` GET on the
   `/export?format=csv` URL (published/link-shared sheet, no auth). Fails loudly on non-200.
2. `train_hallucination_checker(csv, out_dir)` — **TEMPLATE**: if the CSV has `product`/
   `reactants`/`label`, run `ReactionDataLoader.create_dataset` → `stratified_split` →
   `HallucinationClassifier.fit` → `.evaluate` → `.save`; return the saved path. If columns
   are absent, log "template — supply a labeled CSV" and skip (batch still runs with the
   heuristic hallucination mode). Clearly delimited for the other contributor.
3. `read_molecules(txt)` — one SMILES per line; blank/`#`-comment lines skipped.
4. Per molecule → **multiple routes**: take the top-K first-step candidate pathways the
   LLM/agent proposes, solve each to completion, and write each as
   `out/<timestamp>/<safe_molecule>/pathway_<i>.json`. `<timestamp>` is computed once at
   startup (`YYYY-MM-DD_HH-MM-SS`); `<safe_molecule>` is a filesystem-safe slug of the SMILES.
   Per-molecule failures write `error.json` and never abort the batch.

**Multiple-pathways note:** the current `AutoSolver.solve` returns only the first fully-solved
route (DFS). To satisfy the `pathway.json (multiple)` folder shape, the script drives the
solver to emit the top-K candidate routes per molecule (via the first-level candidate
pathways). This is additive; `solve`'s single-route contract is unchanged.

## 7. Error handling

- LLM/tool/sandbox failures degrade gracefully: status codes, `unsolved_leaf`, structured
  sandbox errors, capped agent iterations. No exception escapes a single molecule's solve.
- Batch: per-molecule `try/except` → `error.json`; the run continues.
- Sandbox: wall-clock timeout kills the child; OOM/rlimit → structured error, not a crash.
- Config: missing `--sheet-url`/`--molecules` → argparse error; missing API keys → clear message.

## 8. Testing

- `agents/tools.py`: each executor (valid/invalid SMILES, stability assessment,
  hallucination template stub + injected-classifier verdict, `run_python` happy/error);
  dispatcher unknown-tool path.
- `agents/sandbox.py`: happy path, timeout kill, memory-cap trip, stderr capture,
  temp-dir cleanup.
- `agents/loop.py`: fake `llm_runner` scripting `tool_use` → final answer; asserts tools are
  invoked, iteration cap respected, output shape matches `llm_pipeline`.
- `algorithms/autosolve.py`: `reasoning_mode` wiring with injected `az_runner`/`llm_runner`
  fakes (no network); recursion terminates on AZ-solve; `canonicalize` cycle detection.
- `scripts/run_batch.py`: temp dirs + fake solver → verifies folder structure, top-K files,
  `error.json` on failure; `download_sheet_csv` with a stubbed HTTP response.
- `utils_molecule.canonicalize`: canonical-form + invalid input.

Style: dependency injection + real pharmaceutical molecules (aspirin/paracetamol), **no
monkeypatch/MagicMock** (matches the existing suite). Live-LLM tests marked `slow`/skippable.
`uv run pytest`; `ruff check`/`ruff format`; `ty`; numpydoc docstrings with usage examples.

## 9. Non-goals / templates

- **Hallucination trained-model wiring** — `check_hallucination` tool + `train_hallucination_checker`
  are templates; another contributor completes them. The heuristic path remains fully functional.
- **Google auth** — none; public CSV export URL only.
- **Top-level orchestrator** — scaffolded + flag-wired, minimal impl; `single_step_agent` is
  the fully-built agent mode.
- **Podman/remote sandbox** — protocol seam only; subprocess backend ships.
- No changes to the Flask API / `src/` tree; no changes to the default `pipeline` behavior.

## 10. Required inputs still needed from the user

- The Google Sheets `/export?format=csv` URL.
- The path to the molecules `.txt` file (and confirmation of one-SMILES-per-line format).

Both are wired as required CLI args with placeholders + TODO so the script is runnable the
moment they are supplied; their absence does not block implementation or tests.

## 11. Process

Spec → **4-persona review** (Codex, adversarial reviewer, senior open-source dev, senior
engineer) → revise → implement (TDD) → PR to `dev` (no Claude co-author).
