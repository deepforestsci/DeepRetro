# Design: Unblock CI — Anthropic sampling parameters and langfuse pin

Date: 2026-07-20
Status: Approved (design), not yet planned or implemented
Scope: `src/utils/llm.py`, `tests/requirements_tests.txt`, `tests/variables_test.py`, `tests/test_llm.py`

## Problem

Every open PR fails the same CI job — `test (ubuntu-latest, 3.9)`, the legacy
`tests/` suite. Verified failing on #256, #257, #251, and #240. This is not
PR-specific breakage; it is a main-line environment failure that every branch
inherits.

Two failures are stacked, and the second one hides the first.

### Primary failure (the blocker)

```
FAILED test_adv_prompt.py::test_claude_adv_success - assert 400 == 200
FAILED test_llm.py::test_call_llm_success          - assert 400 == 200
```

Root cause, from the CI log:

```
litellm.exceptions.BadRequestError: AnthropicException -
{"type":"error","error":{"type":"invalid_request_error",
 "message":"`temperature` is deprecated for this model."}}
```

`src/utils/llm.py:150-157` builds `params` with `temperature`, `top_p`, and
`seed` unconditionally. Both failing tests resolve to
`claude-opus-4-20250514` (`tests/variables_test.py:5,7`), which now rejects
`temperature`. The retry at `src/utils/llm.py:190` re-sends identical params,
so it fails the same way twice, and `call_LLM` returns `400, ""` — which is
what the tests assert against.

205 tests pass. Only the two live-API tests fail.

### Secondary failure (noise, not the blocker)

```
litellm/integrations/langfuse/langfuse_prompt_management.py:106
    client = Langfuse(**parameters)
TypeError: __init__() got an unexpected keyword argument 'sdk_integration'
```

`tests/requirements_tests.txt` pins neither `litellm` nor `langfuse`, so CI
resolves both to latest on every run. Langfuse v3 dropped `sdk_integration`
from `Langfuse.__init__`; the resolved litellm still passes it.

This fires **inside the exception handler** for the 400 above — it is
downstream of the real error, not the cause. Pinning langfuse alone will not
turn CI green.

### Note on prior work

`02_Current_State.md` in the SecondBrain vault records this fix as already
made (omit `top_p`, clamp temperature, omit sampling params for
`claude-opus-4-*`, plus regressions named
`test_call_llm_current_opus_model_omits_sampling_parameters`).
`git log --all -S` finds no trace on any branch, and neither the code nor
those tests exist in the tree. That work appears to have been lost from a
worktree without ever being committed. The vault is ahead of git history here.

## Goal

Get CI green with the smallest correct change. Explicitly **not** in scope:
removing live LLM calls from the `tests/` suite. See Known Limitations.

## API facts this design rests on

Verified against the Claude API reference, not assumed:

- `temperature`, `top_p`, and `top_k` are **removed** on Opus 4.7 and Opus 4.8
  — sending any of them returns a 400. The same applies to Sonnet 5 and
  Fable 5.
- Opus 4.6 and Sonnet 4.6 still **accept** sampling parameters.
- `seed` is not a Messages API parameter at all. LiteLLM either drops it or
  passes it through; removing it for Anthropic models is strictly correct.
- `claude-opus-4-20250514` (Claude Opus 4, May 2025) is deprecated with a
  published retirement of 2026-06-15 — already past as of this writing. It is
  still being served (CI receives a 400 parameter error, not a 404), but it is
  living on borrowed time.

Because Opus 4.6/Sonnet 4.6 still accept these parameters, and because
`call_LLM` also routes DeepSeek, OpenAI, Together, and Fireworks models
through the same code path, a blanket removal would be wrong. The gate must
be model-aware.

## Design

### 1. Model-aware sampling-parameter gate — `src/utils/llm.py`

Add a helper alongside the existing `"3-7" in LLM` branch:

```python
# Anthropic removed temperature/top_p/top_k on Opus 4.7+ (400 on send);
# claude-opus-4-* deprecates temperature. `seed` is not an Anthropic param.
_NO_SAMPLING_PREFIXES = ("claude-opus-4", "claude-sonnet-5", "claude-fable-5")


def _accepts_sampling_params(model: str) -> bool:
    return not model.startswith(_NO_SAMPLING_PREFIXES)
```

Applied after `params` is constructed and **after** the existing `"3-7"`
branch, so the thinking-model path keeps its explicit `temperature = 1`:

```python
if not _accepts_sampling_params(LLM):
    params.pop("temperature", None)
    params.pop("top_p", None)
    params.pop("seed", None)
```

Non-Anthropic models are untouched.

Two trade-offs accepted:

**The `claude-opus-4` prefix deliberately over-matches.** It covers the entire
Opus 4 family, including `claude-opus-4-6` — which, per the API facts above,
still *accepts* sampling parameters. This is intentional: stripping a
parameter from a model that tolerates it loses temperature control but is not
an error, whereas sending one to a model that rejects it is a hard 400. Uniform
behaviour across the Opus 4 family is simpler than tracking which point
releases changed. If per-model temperature control on Opus 4.6 turns out to
matter, narrow the prefix to the specific rejecting versions at that point.

**The list is hardcoded** and needs updating as Anthropic ships models. The
alternative considered was an error-driven retry (strip the named parameter on
a `BadRequestError` and retry once), which is self-healing but burns a failed
API call each time and requires parsing provider error strings. The allowlist
was chosen for being explicit and testable without network access.

### 2. Pin `langfuse` and `litellm` — `tests/requirements_tests.txt`

Change `langfuse` to `langfuse<3`, matching the pin already present in
`deepretro/pyproject.toml`. Pin `litellm==1.78.0`, matching `environment.yml:16`
so the conda environment and CI resolve the same version.

This does not fix the failure. It stops the `sdk_integration` TypeError from
burying the real error in the next failure, and stops CI from resolving a
different litellm/langfuse pair on every run — which is how this class of
break recurs.

Related: PR #240 (dependabot, `Update langfuse requirement from <3 to <4 in
/deepretro`) would remove the one pin currently protecting the package path.
It should not be merged as-is.

### 3. Bump the test model — `tests/variables_test.py`

Change `CLAUDE_MODEL` and `CLAUDE_ADV_MODEL` from `claude-opus-4-20250514` to
`claude-opus-4-8`. The gate in (1) already covers `claude-opus-4-8` via the
`claude-opus-4` prefix, so no additional code change is needed.

Rationale: the current model is past its published retirement date. When it
starts returning 404 instead of a parameter error, these tests fail again and
no parameter fix will help.

Correction found during implementation: bumping `variables_test.py` covers
`test_claude_adv_success`, which reads `CLAUDE_ADV_MODEL`, but **not**
`test_call_llm_success` — that test called `call_LLM` with no `LLM` argument and
so inherited the production default in `src/utils/llm.py`, which is still
`claude-opus-4-20250514`. It is now passed `CLAUDE_MODEL` explicitly, so the
bump actually reaches it. The production defaults across `src/` are left alone;
retiring those is a separate decision with a wider blast radius.

### 4. Non-network regression test — `tests/test_llm.py`

Mock `litellm.completion` and assert on the captured kwargs:

- `call_LLM` with an Opus 4 model omits `temperature`, `top_p`, and `seed`.
- `call_LLM` with a DeepSeek model still includes them.

This is the part that stays green independent of Anthropic, secrets, and the
network. Without it the fix is unverifiable in CI and will silently regress —
which is plausibly how the prior attempt was lost.

## Testing

- New mocked regressions in `tests/test_llm.py` (above) — the durable check.
- `python -m pytest tests/ --basetemp=tmp\pytest` locally, per `CLAUDE.md`.
- CI: the two live tests should return 200 once the parameter gate lands,
  assuming valid `ANTHROPIC_API_KEY` in the `testing` environment and that
  `claude-opus-4-8` is available to that key.

## Known limitations

`test_call_llm_success` and `test_claude_adv_success` still make **live,
billed** Anthropic calls from CI. They remain dependent on network, on secrets
being valid, and on Anthropic not deprecating another parameter. This design
does not address that; it was explicitly deferred in favour of a minimal
unblock.

The follow-up — mocking the LLM boundary in `tests/` so CI stops depending on
provider behaviour — is the durable fix. Given that the pinned model is
already past its retirement date, the runway on that follow-up is shorter than
it appears.
