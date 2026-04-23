#!/bin/bash
set -e

REPO_DIR="/mnt/c/Users/gask4/DeepRetro"
MAX_RETRIES=3
LOG_DIR="$REPO_DIR/logs"
VENV_DIR="$REPO_DIR/.testvenv"
CODEX_PROMPT_TEMPLATE="You are a senior Python engineer. PROJECT: DeepRetro FOCUS: Fix failing pytest tests ERRORS: ERRORS_PLACEHOLDER FILE (if available): FILE_PLACEHOLDER CONSTRAINTS: * Do NOT change function signatures * Do NOT modify unrelated logic * Keep code style consistent OUTPUT: Return ONLY a unified diff patch. NO explanations."

cd "$REPO_DIR"

# Ensure logs dir exists
mkdir -p "$LOG_DIR"

# Set up test virtual environment if needed
setup_venv() {
    if [ ! -d "$VENV_DIR" ] || [ ! -f "$VENV_DIR/bin/python" ]; then
        echo "[SETUP] Creating test virtual environment..."
        uv venv "$VENV_DIR" 2>&1
        uv pip install --python "$VENV_DIR/bin/python" pytest pytest-cov 2>&1 | tail -5
        echo "[OK] Virtual environment ready"
    fi
}

PYTEST="$VENV_DIR/bin/python -m pytest"

# Git checkpoint if uncommitted changes
if ! git diff --quiet || [ -n "$(git status --porcelain)" ]; then
    echo "[INFO] Uncommitted changes detected. Creating checkpoint commit..."
    if git add -A && git commit -m "Checkpoint before auto-fix run $(date +'%Y-%m-%d %H:%M:%S')" 2>&1; then
        echo "[OK] Checkpoint created"
    else
        echo "[WARN] Could not create checkpoint, proceeding anyway..."
    fi
fi

run_pytest() {
    $PYTEST -x --tb=short 2>&1
}

echo "=== DeepRetro Auto-Fix Pipeline ==="
echo "Repo: $REPO_DIR"
echo "Max retries: $MAX_RETRIES"
echo ""

setup_venv

retry=0
while [ $retry -lt $MAX_RETRIES ]; do
    echo ""
    echo "--- Attempt $((retry+1))/$MAX_RETRIES ---"
    
    # Run pytest and capture output
    echo "[RUN] pytest..."
    if pytest_output=$(run_pytest 2>&1); then
        echo "[PASS] All tests passed!"
        echo "$pytest_output"
        exit 0
    fi
    
    # Save errors
    error_file="$LOG_DIR/errors_$retry.txt"
    echo "$pytest_output" > "$error_file"
    echo "[FAIL] Tests failed. Saved to $error_file"
    
    # Extract key error lines for context
    echo ""
    echo "--- Failure summary ---"
    tail -20 "$error_file" | grep -E "(FAILED|ERROR|AssertionError|ImportError|NameError)" | head -10
    echo ""
    
    # Build file context
    file_context=""
    if [ -f "$REPO_DIR/deepretro/utils/parse.py" ]; then
        file_context=$(cat "$REPO_DIR/deepretro/utils/parse.py" 2>/dev/null || echo "FILE NOT FOUND")
    else
        file_context="FILE NOT FOUND"
    fi
    
    # Escape special chars for prompt (handle newlines and special chars)
    errors_escaped=$(python3 -c "import sys,json; print(json.dumps(sys.stdin.read()))" < "$error_file" 2>/dev/null | tr -d '"')
    file_escaped=$(python3 -c "import sys,json; print(json.dumps(sys.stdin.read()))" <<< "$file_context" 2>/dev/null | tr -d '"')
    
    # Build prompt
    prompt="${CODEX_PROMPT_TEMPLATE/ERRORS_PLACEHOLDER/$errors_escaped}"
    prompt="${prompt/FILE_PLACEHOLDER/$file_escaped}"
    
    # Call codex CLI with PTY for interactive input
    echo "[CALL] Codex CLI..."
    fix_file="$LOG_DIR/fix_$retry.diff"
    
    codex --yes "$prompt" > "$fix_file" 2>&1
    codex_exit=$?
    
    if [ $codex_exit -ne 0 ]; then
        echo "[WARN] Codex CLI returned exit code $codex_exit"
    fi
    
    echo "[INFO] Codex output saved to $fix_file"
    head -10 "$fix_file"
    
    # Check if codex returned a patch
    if [ ! -s "$fix_file" ]; then
        echo "[FAIL] Codex returned empty output. Giving up."
        exit 1
    fi
    
    # Validate unified diff format
    if ! grep -q "^--- " "$fix_file" || ! grep -q "^+++" "$fix_file"; then
        echo "[WARN] Output doesn't look like a unified diff. Extracting diff portion..."
        python3 -c "
import sys, re
content = open('$fix_file').read()
# Try to find a diff pattern
match = re.search(r'(--- .*\n\+\+\+ .*\n(?:@@ .*\n(?:[+@-].*\n)*)+)', content, re.MULTILINE)
if match:
    print(match.group(1), file=open('$fix_file', 'w'))
else:
    print('[ERROR] Could not extract diff', file=open('$fix_file', 'w'))
" 2>&1
    fi
    
    # Apply patch
    echo "[APPLY] git apply $fix_file..."
    if git apply "$fix_file" 2>&1; then
        echo "[OK] Patch applied successfully"
    else
        echo "[WARN] Patch failed. Retry with --3way..."
        if git apply --3way "$fix_file" 2>&1; then
            echo "[OK] Patch applied via 3-way merge"
        else
            echo "[FAIL] Patch could not be applied. Saving and retrying..."
            cp "$fix_file" "$LOG_DIR/failed_fix_$retry.diff"
            retry=$((retry+1))
            continue
        fi
    fi
    
    # Stage the changes
    git add -A
    
    retry=$((retry+1))
done

echo ""
echo "=== Exhausted $MAX_RETRIES retries without full pass ==="
echo "Last error file: $error_file"
echo "Last fix attempt: $fix_file"
exit 1