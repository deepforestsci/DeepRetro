# Dev versus Main Functionality Report Design

## Objective

Create a comprehensive Markdown report explaining every net functionality on
the upstream GitHub `dev` branch that is not present on `main`, using important
code snippets and complete file-level traceability. Open the finished report in
the user's default desktop application for immediate review.

## Authoritative Comparison

The report is based on freshly verified upstream GitHub branch heads:

- repository: `https://github.com/deepforestsci/DeepRetro`
- `main`: `dc8fcad1585650c58d4e4adafca710ca007fb509`
- `dev`: `0575a1ffaa2c18f3681744cf1494c285090a257b`
- merge base: `dc8fcad1585650c58d4e4adafca710ca007fb509`

The effective comparison is `upstream/main...upstream/dev`. It contains 47
commits, 51 changed files, 18,559 insertions, and 97 deletions. Local uncommitted
changes are excluded.

## Report Approach

Use a hybrid structure:

1. an executive summary of the change surface;
2. a system map preserving the `src/` runtime versus `deepretro/` package split;
3. thematic functionality sections with important, short code snippets;
4. operational implications and known limitations;
5. a complete file-by-file traceability appendix covering all 51 paths; and
6. a chronological appendix covering all 47 commits, including merges and
   reversions.

The report explains the final net behavior. Intermediate code that was later
reverted is identified in the chronology but is not presented as current
functionality.

## Functionality Groups

The report will organize the diff into these groups:

- agentic single-step retrosynthesis and tool execution;
- hardened local Python sandboxing;
- package-level autosolve orchestration and recursive route handling;
- AiZynthFinder caching and agent tools;
- LLM provider selection, adaptive thinking, prompt caching, timeouts, and
  retry controls;
- hallucination checker adaptation, ML training, and dataset additions;
- batch solving and command-line entry points;
- SA/SCScore pathway scoring and report generation;
- reaction validation, featurization, and molecule utility corrections;
- documentation and package dependency updates;
- CI workflow targeting for `dev`; and
- runtime compatibility changes under `src/`.

## Code Snippet Rules

Snippets will be copied from the remote-ref diff, not the dirty checkout. Each
snippet will be short enough to explain one contract or control flow and will
name its source path. Long implementations, generated lockfiles, notebooks, and
CSV contents will be summarized rather than reproduced.

## Completeness Checks

Before delivery:

- compare the appendix path list with `git diff --name-status` and require all
  51 changed paths to appear;
- compare the chronology with `git log upstream/main..upstream/dev` and require
  all 47 commits to appear;
- verify every thematic claim against the final remote-ref source or diff;
- scan for placeholders and ambiguous language;
- run `git diff --check` for the report; and
- confirm the opened path is the final Markdown file.

## Output

Write the report to:

`docs/dev-vs-main-functionality-report.md`

After verification, open it with the Windows default application via
PowerShell's `Start-Process`, the local equivalent of a shell `open` command.
