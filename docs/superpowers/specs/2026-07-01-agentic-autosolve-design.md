# Agentic tool-calling retrosynthesis in `deepretro/`

**Status:** Approved design, revised after 4-persona review (Codex + adversarial + senior OSS
maintainer + senior engineer). Ready to implement.
**Branch:** `feat/agentic-autosolve` (off `dev`) → single PR to `dev` (no Claude co-author).
**Date:** 2026-07-01

## 1. Goal

Adapt DeepRetro's recursive retrosynthesis algorithm to run as an **agentic loop with tool
calling**, inside the `deepretro/` package. Three deliverables:

1. Port the recursive loop + metadata into `deepretro/` (built on the unmerged
   `origin/autosolve-pipeline` branch, rewired to current `dev` APIs).
2. An agentic layer where the LLM can call tools — **molecule validity**, **stability**,
   **hallucination** (template), and **write-and-run-Python-in-a-sandbox** — with structured
   tool calling and sandbox code-execution both selectable.
3. A batch script: download a CSV from Google Sheets → train the hallucination checker
   (template) → run molecules from a `.txt` file → dump `folder/<timestamp>/<molecule>/pathway_<i>.json`.

## 2. Current state (verified on `dev`)

- **LLM layer** (100% `litellm`): `utils/llm_interface.py` (`LLMInterface.call()` returns a
  **text-only** `LLMResponse(status_code, text)`), `utils/llm.py` (`call_LLM`, `llm_pipeline`
  — single-step retrosynthesis: rising-temperature retry loop → parse → `validate_split_json`
  → `validity_check` → optional `filter_stable_pathways` → optional `filter_hallucination_pathways`),
  `utils/llm_helpers.py` (`build_completion_params` — the one place assembling
  `litellm.completion(**params)`; **has no `tools`/`tool_choice` path today**).
- **Checks**: `algorithms/stability_checker.py::check_molecule_stability`,
  `algorithms/hallucination_checker.py::calculate_hallucination_score` (heuristic),
  `utils/utils_molecule.py::is_valid_smiles` / `validity_check`,
  `algorithms/pipeline_checks.py::hallucination_checker` (order-preserving pathway filter).
- **ML hallucination**: `models/hallucination_classifier.py::HallucinationClassifier`
  (XGBoost via DeepChem `GBDTModel`; `fit`, `evaluate`, `predict_probability`,
  `save(save_dir)`, `load(save_dir)`; exposes `.featurizer`, `.threshold`),
  `models/hallucination_helpers.py::MLChecker` + `resolve_hallucination(mode, classifier)`
  (returns a `Callable[[str, list], tuple[int, list]]` or `None`).
- **Metadata**: `metadata.py::recommend_reaction_metadata`. **AZ**: `utils/az.py::run_az`
  (can raise `FileNotFoundError`/`ImportError` before returning a status tuple).
  **Parse**: `utils/parse.py::format_output`. **Training path**: `data/loader.py::ReactionDataLoader`
  (cols `product`/`reactants`/`label`) + `stratified_split`, `featurizers/reactionstep.py`.
- `requests` is a direct dependency (via `aizynthfinder`/tooling); `httpx` is only transitive.

Missing on `dev`: the recursive loop, any tool-calling/agent layer, any sandbox. A dead `TOOLS`
list sits unused in `utils/variables.py`; `origin/agentic-solution` has a one-tool prototype
(`agentic_call_LLM` + `src/utils/tools.py`) — a **pattern reference only**; 3 of the 4 tools and
the entire sandbox here are net-new.

## 3. Reuse + adaptation layer

`origin/autosolve-pipeline:deepretro/algorithms/autosolve.py` (~636 lines) implements `AutoSolver`
(recursive `solve`, `single_step`, `parse`, `add_metadata`, `autosolve`, `run_llm`) with DI, plus
~660 lines of tests. It targets the older `autosolve-foundations` API and is **not** cleanly
rebasable. Reuse it as the structural basis with this bounded adaptation on `dev`:

| Branch symbol | On `dev` | Adaptation |
|---|---|---|
| `hallucination_helpers.filter_with_checker(molecule, pathways, explanations, confidence, checker)` | absent | **Add it** to `dev`'s `hallucination_helpers.py`: `checker=None → passthrough`; else `status, kept = checker(molecule, pathways)`; `status != 200 → ([],[],[])`; realign `explanations`/`confidence` to `kept` by **order-preserving index walk** (dev checkers preserve order). Unit-tested. |
| `hallucination_helpers.HallucinationPredictor` (protocol w/ `predict_single`) | absent | Drop it. Type `hallucination_classifier: Any` (dev contract: an object with `predict_probability`/`threshold`/`featurizer`, or a `str`/`Path` — whatever `resolve_hallucination("ml", …)` accepts). |
| `utils.typing.HallucinationChecker` | absent | Add alias `HallucinationChecker = Callable[[str, list], tuple[int, list]]`. |
| `utils_molecule.canonicalize(smiles) -> str` | absent | Port the trivial RDKit wrapper (`MolToSmiles(MolFromSmiles(s), canonical=True)`; `""` on invalid). Unit-tested. |
| `run_az`, `llm_pipeline`, `format_output`, `recommend_reaction_metadata`, recommender types | present | no change |

## 4. Architecture

```
deepretro/
├── algorithms/autosolve.py     ported AutoSolver (rewired) + solve_mode/tool_backend + solve_multiple
├── agents/                     NEW subpackage
│   ├── __init__.py
│   ├── tools.py                tool registry (OpenAI function schemas) + executors + dispatcher
│   ├── sandbox.py              Sandbox protocol + SubprocessSandbox (secret-scrubbed, rlimited)
│   └── loop.py                 agentic_single_step() (full) + agentic_orchestrator() (NotImplementedError)
├── utils/utils_molecule.py     + canonicalize()
├── utils/typing.py             + HallucinationChecker alias
├── utils/llm_helpers.py        build_completion_params gains optional tools/tool_choice
├── utils/llm.py                extract validity+stability filters as reusable helpers (already ~separate)
└── models/hallucination_helpers.py  + filter_with_checker() (with realignment)
scripts/run_batch.py            NEW batch driver
```

### 4.1 `AutoSolver` flags and post-filtering (item: "both, flag-selectable")

Constructor gains:
- `solve_mode: {"pipeline", "single_step_agent", "orchestrator"} = "pipeline"` (renamed from
  `reasoning_mode` to avoid colliding with `enable_thinking`/reasoning-effort concepts).
- `tool_backend: {"structured", "sandbox"} = "structured"` (renamed value `tools`→`structured`;
  **ignored in `pipeline` mode**, documented as such).
- `sandbox: Sandbox | None = None` (injectable; defaults to `SubprocessSandbox`).

`run_llm(molecule)` is the **single owner of post-filtering** for both modes, so the agent path
can never bypass the safety filters:

- **`pipeline`** (default, unchanged behavior): call `llm_pipeline(..., stability_check=self.stability_check,
  hallucination_check=False)` (validity + stability owned by `llm_pipeline`; hallucination turned
  **off** there to avoid a double pass), then apply the resolved hallucination checker via
  `filter_with_checker`.
- **`single_step_agent`** (fully built): call `agents.loop.agentic_single_step(...)` to get **raw**
  parsed `(pathways, explanations, confidence)`, then apply the shared safety filter
  `_apply_safety_filters` = `validity_check` → `filter_stable_pathways` (if `stability_check`) →
  `filter_with_checker` (hallucination). Same output contract as `llm_pipeline`.
- **`orchestrator`**: `agents.loop.agentic_orchestrator` **raises `NotImplementedError`** with a
  clear message. Interface + flag wiring exist; not offered in the batch CLI. (Non-goal §9.)

Hallucination is owned in exactly one place (`filter_with_checker` in `run_llm`) for both modes.

### 4.2 Tools (`agents/tools.py`) — OpenAI function-call format

litellm's tool calling normalizes to OpenAI format, so schemas are stored as
`{"type": "function", "function": {"name", "description", "parameters": <json-schema>}}` (not
Anthropic-native `input_schema`). Registry + `execute_tool(name, input) -> dict` dispatcher +
`get_tool_schemas(tool_backend) -> list`.

| tool | signature | wraps | returns |
|---|---|---|---|
| `validate_smiles` | `(smiles, context?)` | `is_valid_smiles` + `canonicalize` | `{valid, canonical_smiles, error, context}` |
| `check_stability` | `(smiles)` | `check_molecule_stability` | `{assessment, stability_score, issues}` |
| `check_hallucination` | `(product, reactants)` | the checker from `resolve_hallucination(mode, clf)` | verdict from the resolved checker. **TEMPLATE seam:** the executor is bound to the resolved checker at registry-build time; in `heuristic` mode it uses the working heuristic checker (never a no-op), in `ml` mode it uses the injected classifier, in `none` mode the tool is **omitted from the registry**. The ML wiring is what another contributor completes. |
| `run_python` | `(code)` | `agents/sandbox.py` | `{ok, stdout, stderr, error}` — RDKit + `deepretro.algorithms` importable |

`tool_backend="structured"` exposes `validate_smiles`/`check_stability`/`check_hallucination`;
`tool_backend="sandbox"` exposes `run_python` (+ the read-only checks). The registry is built per
`AutoSolver`, binding `check_hallucination` to the resolved checker so the ML classifier swaps in
without touching tool code.

### 4.3 Sandbox (`agents/sandbox.py`) — hardened subprocess (the chosen backend)

```python
class Sandbox(Protocol):
    def run(self, code: str, *, timeout_s: float = 10.0, mem_mb: int = 2048) -> SandboxResult: ...

class SubprocessSandbox:   # default
    # - writes code to a temp file in a fresh temp cwd (cleaned up after)
    # - MANDATORY secret scrubbing: child runs with a minimal allowlist env
    #   (PATH, LANG, and nothing else) — NO ANTHROPIC_API_KEY/OPENAI_API_KEY/etc. inherited
    # - POSIX rlimits (best-effort, per-limit try/except): RLIMIT_AS (mem_mb*1024*1024, with a
    #   floor of 1024 MB so RDKit/deepretro imports don't MemoryError), RLIMIT_CPU (ceil timeout),
    #   RLIMIT_NPROC, RLIMIT_FSIZE
    # - hard wall-clock timeout: kill the child (and its process group) on expiry
    # - best-effort network isolation: wrap with `unshare -rn` when available (Linux); on macOS
    #   there is NO network isolation — documented ceiling
    # - captured stdout/stderr; returns SandboxResult(ok, stdout, stderr, error)
```

`ponytail:` comment names the ceiling: secret-scrubbing + rlimits + best-effort `unshare -n` is
adequate for **LLM-generated, semi-trusted** chem code, but is **not** a true jail; upgrade path is
a `PodmanSandbox` (`--network=none`, read-only rootfs) implementing the same `Sandbox` protocol.
**Threat-model statement in the module + docs:** untrusted Google-Sheet CSV / molecule-file text is
NEVER interpolated into `run_python` code; the agent generates code from its own reasoning; the
primary mitigations are secret-scrubbing (no key exfiltration) + resource caps + (optional) network
isolation. `RLIMIT_AS` is unenforced on macOS — treated as best-effort there.

### 4.4 Agent loop (`agents/loop.py`) — widened tool-response contract

The text-only `LLMInterface.call()` cannot carry `tool_calls`, so the loop calls
`litellm.completion` **directly** through an extended `build_completion_params(..., tools=…,
tool_choice=…)` and reads the **full** message (`content` + `tool_calls`) and `finish_reason`:

```python
def agentic_single_step(
    molecule, model, *, tool_registry, tool_backend, max_iterations=6,
    llm_runner=None, enable_thinking=False, max_output_tokens=None,
) -> tuple[list[Pathway], list[str], list[float]]:
    # system prompt = the pipeline's retrosynthesis prompt + "you may call tools to check
    #   candidates; when done, emit the SAME <json> payload the pipeline expects"
    # loop while message.tool_calls is non-empty and iter < max_iterations:
    #     append assistant turn (with tool_calls); for each tool_call:
    #         execute_tool(name, json.loads(arguments)) -> append role="tool" result msg
    # on a final (no-tool_calls) message: reuse parse_response/validate_split_json ->
    #     (pathways, explanations, confidence)   [same shape as llm_pipeline]
    # if max_iterations hit with no final answer: return ([], [], []) and log the reason
```

Details fixed per review:
- Loop condition is **non-empty `message.tool_calls`** (litellm normalizes Anthropic to
  `finish_reason=="tool_calls"`); no reliance on a literal `"tool_use"`.
- `enable_thinking=False` by default in the tool loop — extended thinking + multi-turn tool use
  returns Anthropic HTTP 400 on the second turn unless the signed thinking block is preserved.
- `max_iterations` exhausted → `([], [], [])` → `run_llm` → `unsolved_leaf`; asserted in a
  fake-runner iteration-cap test.
- `llm_runner` injectable: a fake scripting `tool_calls` then a final `<json>` answer — no network.

## 5. Data flow (recursion terminates when AZ solves)

```
autosolve(smiles):  route_tree, solved = solve(smiles); parse; add_metadata

solve(molecule, visited=None, depth=0):
    m = _clean_smiles(molecule);  if not m: return unsolved_leaf, False
    canonical = canonicalize(m)
    if depth >= max_depth: return unsolved_leaf, False          # safety guard
    branch_visited = set(visited or ())                         # PER-BRANCH COPY (path-only dedup)
    if canonical in branch_visited: return unsolved_leaf, False # cycle guard
    branch_visited.add(canonical)
    try: az_solved, az_routes = az_runner(m, az_model)          # wrapped: FileNotFoundError/ImportError -> unsolved_leaf
    except (FileNotFoundError, ImportError, ...): return unsolved_leaf, False
    if az_solved and az_routes: return az_routes[0], True       # ✅ AZ found building blocks -> branch ends
    pathways, expl, conf = run_llm(m)                           # pipeline OR agent, per solve_mode; single-owner filtering
    for i, pathway in enumerate(pathways, best-first):
        children, all_solved = _solve_pathway(pathway, branch_visited, depth)   # recurse each reactant
        if all_solved: return reaction_tree(m, children, [conf[i]]), True
    return best_effort_partial, False
```

A branch stops **only** when AZ solves the node. A pathway is solved iff every leaf bottoms out at
an AZ-solved node. `visited` is a **per-branch copy** (path-only dedup — the `autosolve-pipeline`
semantics, correct for a tree; this differs from `rec_prithvi`'s global set, and we adopt the
per-branch copy deliberately). `max_depth`/cycle guards prevent infinite recursion.

### 5.1 Top-K routes (`solve_multiple`)

`solve()` returns a single route. Add:
```python
def solve_multiple(self, smiles, k=3) -> list[tuple[dict, bool]]:
    # AZ first: if AZ solves, return [(az_route, True)]
    # else: pathways,_,conf = run_llm(smiles); for each of the top-k pathways, solve each reactant
    #       with its OWN fresh visited set -> reaction_tree; return up to k (route, solved) tuples
```
Each candidate route is independent (fresh visited), producing the `pathway_<i>.json` deliverable.

## 6. Batch script (`scripts/run_batch.py`)

```
python scripts/run_batch.py \
  --sheet-url <google-sheets /export?format=csv URL>   # REQUIRED (public export, no auth)
  --molecules <path/to/molecules.txt>                  # REQUIRED (one SMILES per line)
  --out <output-dir> --solve-mode {pipeline,single_step_agent} \  # orchestrator NOT offered
  --tool-backend {structured,sandbox} --model <litellm id> \
  --az-model <default USPTO> --classifier <optional trained-model dir> --top-k <int, default 3>
```

Steps:
1. `download_sheet_csv(export_url, dest)` — `requests.get` on the `/export?format=csv` URL
   (published/link-shared, no auth); raise on non-200.
2. `train_hallucination_checker(csv, save_dir) -> str | None` — **TEMPLATE**: if the CSV has
   `product`/`reactants`/`label`, `ReactionDataLoader.create_dataset` → `stratified_split` →
   `clf = HallucinationClassifier(model_dir=save_dir)` → `clf.fit(train)` → `clf.evaluate(test)`
   → `clf.save(save_dir)` (explicit `save_dir`; the function returns `save_dir`). If columns are
   absent: log "template — supply a labeled CSV", return `None`, batch continues in `heuristic`
   hallucination mode. Clearly delimited for the other contributor.
3. `read_molecules(txt)` — one SMILES per line; blank/`#` lines skipped.
4. Per molecule → `AutoSolver(...).solve_multiple(smiles, k=top_k)`; parse + add_metadata each
   route; write `out/<timestamp>/<safe_molecule>/pathway_<i>.json`. `<timestamp>` computed once
   at startup (`YYYY-MM-DD_HH-MM-SS`); `<safe_molecule>` = filesystem-safe slug. Per-molecule
   `try/except` → `error.json`; the batch never aborts.

## 7. Error handling

- `solve` wraps `az_runner`/AZ internals (`FileNotFoundError`/`ImportError`) → `unsolved_leaf`,
  so no exception escapes a single molecule's solve. LLM/tool/sandbox failures degrade to status
  codes / `unsolved_leaf` / structured errors. Agent iterations capped.
- Batch: per-molecule `try/except` → `error.json` backstop; run continues.
- Sandbox: wall-clock timeout kills the child group; OOM/rlimit → structured error.
- Missing `--sheet-url`/`--molecules` → argparse error; missing API keys → clear message.

## 8. Testing

- `models/hallucination_helpers.filter_with_checker`: realignment correctness (kept subset keeps
  matching expl/conf; `status!=200 → empty`; `checker=None → passthrough`).
- `utils_molecule.canonicalize`: canonical form + invalid input.
- `agents/tools.py`: each executor (valid/invalid SMILES; stability assessment; `check_hallucination`
  heuristic verdict + injected-classifier verdict + omitted in `none` mode; `run_python` happy/error);
  unknown-tool dispatch; schemas are OpenAI function format.
- `agents/sandbox.py`: happy path; timeout kill; stderr capture; temp-dir cleanup; **secret scrub**
  (assert `ANTHROPIC_API_KEY` set in parent is NOT visible to child); memory-cap test marked
  **Linux-only** (`RLIMIT_AS` unenforced on macOS).
- `agents/loop.py`: fake `llm_runner` scripting `tool_calls` → final answer; asserts tools invoked,
  output shape matches `llm_pipeline`, iteration cap → `([],[],[])`.
- `algorithms/autosolve.py`: `solve_mode` wiring with injected `az_runner`/`llm_runner` fakes (no
  network); AZ-solve terminates recursion; per-branch visited/cycle detection; `run_llm` applies
  safety filters in agent mode; `solve_multiple` yields ≤k independent routes; AZ raising →
  `unsolved_leaf`; `orchestrator` → `NotImplementedError`.
- `scripts/run_batch.py`: temp dirs + fake solver → folder structure, top-K files, `error.json`;
  `download_sheet_csv` with a stubbed `requests` response.

Style: dependency injection + real pharmaceutical molecules (aspirin/paracetamol), **no
monkeypatch/MagicMock**; live-LLM tests marked `slow`/skippable. `uv run pytest`; `ruff
check`/`ruff format`; `ty`; numpydoc docstrings with usage examples.

## 9. Non-goals / templates

- **Hallucination trained-model wiring** — `check_hallucination` (ml mode) + `train_hallucination_checker`
  are templates for another contributor. The heuristic path is fully functional throughout.
- **Google auth** — none; public CSV export URL only.
- **Top-level orchestrator** — flag + interface only; `agentic_orchestrator` raises
  `NotImplementedError`; not offered in the batch CLI.
- **Podman/remote sandbox** — protocol seam only; hardened subprocess ships.
- No changes to the Flask API / `src/` tree; the default `pipeline` behavior is unchanged.

## 10. Required inputs still needed from the user

- The Google Sheets `/export?format=csv` URL.
- The molecules `.txt` path (one SMILES per line).

Wired as required CLI args with placeholders + TODO; absence blocks neither implementation nor tests.

## 11. Process

Spec (this doc) → 4-persona review **(done; incorporated)** → implement (TDD, logically-ordered
commits) → single PR to `dev` (no Claude co-author). Reviewers suggested stacked PRs; deferring to
the single-PR request, splittable on ask.
```
Implementation order (commits): (1) adaptation layer — canonicalize, HallucinationChecker alias,
filter_with_checker, extracted safety filters + tests; (2) AutoSolver port + solve_mode/solve_multiple
+ tests; (3) agents/ subpackage — tools, sandbox, loop + tests; (4) scripts/run_batch.py + training
template + tests.
```
