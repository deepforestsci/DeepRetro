#!/bin/bash
set -e

REPO_DIR="/mnt/c/Users/gask4/DeepRetro"
MAX_RETRIES=3
LOG_DIR="$REPO_DIR/logs"
CODEX_PROMPT_TEMPLATE="You are a senior Python engineer. PROJECT: DeepRetro FOCUS: Fix failing pytest tests ERRORS: ERRORS_PLACEHOLDER FILE (if available): FILE_PLACEHOLDER CONSTRAINTS: * Do NOT change function signatures * Do NOT modify unrelated logic * Keep code style consistent OUTPUT: Return ONLY a unified diff patch. NO explanations."

cd "$REPO_DIR"

# Ensure logs dir exists
mkdir -p "$LOG_DIR"

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
    pytest -x --tb=short 2>&1
}

echo "=== DeepRetro Auto-Fix Pipeline ==="
echo "Repo: $REPO_DIR"
echo "Max retries: $MAX_RETRIES"
echo ""

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
    tail -20 "$error_file" | grep -E "(FAILED|ERROR|AssertionError|ImportError)" | head -10
    echo ""
    
    # Build file context
    file_context=""
    if [ -f "$REPO_DIR/deepretro/utils/parse.py" ]; then
        file_context=$(cat "$REPO_DIR/deepretro/utils/parse.py" 2>/dev/null || echo "FILE NOT FOUND")
    else
        file_context="FILE NOT FOUND"
    fi
    
    # Escape special chars for prompt
    errors_escaped=$(sed 's/"/\\"/g; s/`/\\`/g' "$error_file")
    file_escaped=$(sed 's/"/\\"/g; s/`/\\`/g' <<< "$file_context")
    
    # Build prompt
    prompt="${CODEX_PROMPT_TEMPLATE/ERRORS_PLACEHOLDER/$errors_escaped}"
    prompt="${prompt/FILE_PLACEHOLDER/$file_escaped}"
    
    # Call codex CLI
    echo "[CALL] Codex CLI..."
    fix_file="$LOG_DIR/fix_$retry.diff"
    
    codex "$prompt" > "$fix_file" 2>&1
    codex_exit=$?
    
    if [ $codex_exit -ne 0 ]; then
        echo "[WARN] Codex CLI returned exit code $codex_exit"
    fi
    
    echo "[INFO] Codex output saved to $fix_file"
    
    # Check if codex returned a patch
    if [ ! -s "$fix_file" ]; then
        echo "[FAIL] Codex returned empty output. Giving up."
        exit 1
    fi
    
    # Enforce unified diff format - if output contains more than diff, extract just the diff
    if grep -q "^--- " "$fix_file" && grep -q "^+++" "$fix_file"; then
        echo "[INFO] Valid diff detected in Codex output"
    else
        echo "[WARN] Output doesn't look like a unified diff. Attempting to extract..."
        # Try to extract diff portion
        sed -n '/^--- /,/^@@ /p' "$fix_file" > "${fix_file}.tmp" 2>/dev/null || true
        if [ -s "${fix_file}.tmp" ]; then
            mv "${fix_file}.tmp" "$fix_file"
        fi
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
            echo "[FAIL] Patch could not be applied. Saving and trying again..."
            cp "$fix_file" "$LOG_DIR/failed_fix_$retry.diff"
            retry=$((retry+1))
            echo "[INFO] Retrying with more context..."
            continue
        fi
    fi
    
    # Verify patch makes sense - stage the changes
    git add -A
    
    retry=$((retry+1))
done

echo ""
echo "=== Exhausted $MAX_RETRIES retries without full pass ==="
echo "Last error file: $error_file"
echo "Last fix attempt: $fix_file"
exit 1