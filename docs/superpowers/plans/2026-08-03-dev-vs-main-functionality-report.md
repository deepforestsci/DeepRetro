# Dev versus Main Functionality Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce and open a comprehensive Markdown explanation of every net functionality added on upstream `dev` relative to upstream `main`.

**Architecture:** Treat remote Git refs as immutable input data and ignore the dirty checkout. Build a thematic explanation from the final diff, then prove completeness with file and commit appendices generated from the same refs.

**Tech Stack:** Git, PowerShell, Markdown, Python AST for read-only symbol inventory, GitHub branch/compare URLs.

---

### Task 1: Pin and inventory the upstream comparison

**Files:**
- Read: remote refs only
- Create: `docs/dev-vs-main-functionality-report.md`

- [ ] **Step 1: Verify remote branch heads**

Run:

```powershell
git fetch upstream main dev
git ls-remote upstream refs/heads/main refs/heads/dev
```

Expected: `main` is `dc8fcad1585650c58d4e4adafca710ca007fb509` and `dev` is `0575a1ffaa2c18f3681744cf1494c285090a257b`.

- [ ] **Step 2: Capture the complete change inventory**

Run:

```powershell
git diff --name-status upstream/main...upstream/dev
git diff --stat upstream/main...upstream/dev
git log --reverse --format='%H`t%ad`t%s' --date=short upstream/main..upstream/dev
```

Expected: 51 changed files and 47 commits.

### Task 2: Derive functionality groups from final source

**Files:**
- Read: all 51 paths from `upstream/main...upstream/dev`
- Create: `docs/dev-vs-main-functionality-report.md`

- [ ] **Step 1: Read final added modules from the remote ref**

Use `git show upstream/dev:<path>` for new source files and `git diff upstream/main...upstream/dev -- <path>` for modified files. Do not read equivalent paths from the working tree.

- [ ] **Step 2: Map each source change to user-visible or developer-visible behavior**

Cover agent tools, sandboxing, autosolve, AiZynthFinder caching, LLM controls, hallucination training, batch execution, scoring/reporting, validation utilities, CI, docs, datasets, scripts, and runtime compatibility.

- [ ] **Step 3: Select representative code snippets**

Copy short snippets directly from `upstream/dev`. Each snippet must name its source path and demonstrate one contract, decision point, or data shape. Summarize generated locks, CSV rows, and notebook cells instead of reproducing them.

### Task 3: Write the report

**Files:**
- Create: `docs/dev-vs-main-functionality-report.md`

- [ ] **Step 1: Write comparison metadata and executive summary**

Include repository links, exact SHAs, merge base, commit/file counts, and diff statistics.

- [ ] **Step 2: Explain each functionality group**

Preserve the distinction between the Flask runtime under `src/` and the refactored package under `deepretro/`. For each group, explain what was added, how it works, why it matters, key files, important snippets, and limitations.

- [ ] **Step 3: Add complete traceability appendices**

List all 51 changed paths with status and functional purpose. List all 47 commits chronologically, noting merge commits and reversions without presenting superseded code as current.

### Task 4: Verify completeness and open the artifact

**Files:**
- Verify: `docs/dev-vs-main-functionality-report.md`

- [ ] **Step 1: Verify path and commit coverage**

Compare the report's file appendix against `git diff --name-status upstream/main...upstream/dev` and its commit appendix against `git log upstream/main..upstream/dev`. Expected: zero missing paths and zero missing commits.

- [ ] **Step 2: Verify document integrity**

Run:

```powershell
Select-String -Path docs\dev-vs-main-functionality-report.md -Pattern '\b(TBD|TODO|FIXME|XXX)\b'
git diff --check -- docs/dev-vs-main-functionality-report.md
```

Expected: no placeholders and no whitespace errors.

- [ ] **Step 3: Open the finished report**

Run:

```powershell
Start-Process -FilePath (Resolve-Path 'docs\dev-vs-main-functionality-report.md')
```

Expected: the Markdown report opens in the Windows default application.
