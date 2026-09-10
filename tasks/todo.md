# Autosolve metadata

- [x] Trace AZ success flags, attribution, parser, and batch paths.
- [x] Mark AZ subtrees in all entry points without mutating runner payloads.
- [x] Preserve LLM/AZ reaction attribution in parsed steps; distinguish AZ-only completion from mixed completion.
- [x] Add mocked regression coverage and document output semantics.
- [x] Run affected tests, lint/type/docs checks, and independent review.
- [x] Commit, push, and open a PR against `dev`.

Acceptance: AZ-only success sets both flags; completed mixed routes set only
`solved`; incomplete routes do not set `az_solved`. Every autosolve reaction
has `solved_by`, and terminal molecules do not become phantom reaction steps.

Validation: final autosolve/parser tests 88 passed; initial autosolve/parser/batch
run 92 passed. Ruff lint/format and Sphinx warnings-as-errors build pass.
Changed source files pass ty. Full-package ty reports the existing optional
AiZynthFinder import assignment at utils/az.py:28; reproduced on HEAD source.
Independent source review found no blocking issues. No live model calls used.

PR: https://github.com/deepforestsci/DeepRetro/pull/281 (base: dev).
