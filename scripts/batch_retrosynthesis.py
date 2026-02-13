"""
Batch retrosynthesis runner - processes SMILES from a file.

Features:
- Automatic resume: tracks progress and continues from last completed on re-run
- Progress tracking: saves state to .progress file

Usage:
    python scripts/batch_retrosynthesis.py
    python scripts/batch_retrosynthesis.py --input data/SMILES.txt
    python scripts/batch_retrosynthesis.py --hallucination-check --stability-check
    python scripts/batch_retrosynthesis.py --reset  # Start fresh, ignore progress
"""

import os
import sys
import argparse
import time
import json
import hashlib
from pathlib import Path
from datetime import datetime

# Add project root to path
import rootutils
root_dir = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

from src.main import main as run_retrosynthesis
from src.cache import clear_cache_for_molecule


PROGRESS_DIR = Path(root_dir) / "data" / ".batch_progress"
RESULTS_DIR = Path(root_dir) / "data" / "batch_results"


def get_progress_file(input_file: str) -> Path:
    """Get progress file path based on input file hash."""
    # Create unique progress file per input file
    file_hash = hashlib.md5(input_file.encode()).hexdigest()[:8]
    PROGRESS_DIR.mkdir(parents=True, exist_ok=True)
    return PROGRESS_DIR / f"progress_{file_hash}.json"


def load_progress(input_file: str) -> dict:
    """Load progress from file."""
    progress_file = get_progress_file(input_file)
    if progress_file.exists():
        with open(progress_file, 'r') as f:
            return json.load(f)
    return {"completed": [], "failed": [], "last_index": -1}


def save_progress(input_file: str, progress: dict):
    """Save progress to file."""
    progress_file = get_progress_file(input_file)
    progress["last_updated"] = datetime.now().isoformat()
    with open(progress_file, 'w') as f:
        json.dump(progress, f, indent=2)


def reset_progress(input_file: str):
    """Reset progress file."""
    progress_file = get_progress_file(input_file)
    if progress_file.exists():
        progress_file.unlink()
        print(f"Progress reset: {progress_file}")


def load_smiles(filepath: str) -> list:
    """Load SMILES from a text file (one per line)."""
    smiles_list = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                smiles_list.append(line)
    return smiles_list


def save_result_json(smiles: str, result: dict, output_dir: Path):
    """Save result JSON to file."""
    # Create safe filename from SMILES (hash if too long)
    if len(smiles) > 50:
        filename = f"{hashlib.md5(smiles.encode()).hexdigest()[:8]}.json"
    else:
        # Replace invalid filename characters
        safe_smiles = smiles.replace("/", "_").replace("\\", "_").replace(":", "_")
        filename = f"{safe_smiles}.json"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / filename
    
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)
    
    return output_path


def run_batch(
    input_file: str = "data/SMILES.txt",
    llm: str = "claude-sonnet-4-20250514",
    az_model: str = "USPTO",
    stability_check: bool = False,
    hallucination_check: bool = True,
    hallucination_method: str = "rule_based",
    use_protecting_groups: bool = False,
    limit: int = None,
    skip: int = 0,
    reset: bool = False,
    clear_cache: bool = False
):
    """Run batch retrosynthesis on SMILES file with auto-resume."""
    
    # Handle reset
    if reset:
        reset_progress(input_file)
    
    # Create results directory with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    results_subdir = RESULTS_DIR / timestamp
    results_subdir.mkdir(parents=True, exist_ok=True)
    
    # Load SMILES
    smiles_list = load_smiles(input_file)
    total_smiles = len(smiles_list)
    print(f"Loaded {total_smiles} SMILES from {input_file}")
    print(f"Results will be saved to: {results_subdir}")
    
    # Clear cache for all molecules if requested
    if clear_cache:
        print(f"\n*** CLEARING CACHE for {total_smiles} molecules ***")
        for smiles in smiles_list:
            clear_cache_for_molecule(smiles)
        print("Cache cleared!\n")
    
    # Load progress
    progress = load_progress(input_file)
    completed_smiles = set(progress.get("completed", []))
    
    # Check for resume
    if completed_smiles:
        print(f"\n*** RESUMING: {len(completed_smiles)} already completed ***")
        print(f"    (use --reset to start fresh)\n")
    
    # Apply skip and limit
    if skip > 0:
        smiles_list = smiles_list[skip:]
        print(f"Skipped first {skip}, remaining: {len(smiles_list)}")
    
    if limit:
        smiles_list = smiles_list[:limit]
        print(f"Limited to first {limit}")
    
    # Filter out already completed
    pending_smiles = [(i, s) for i, s in enumerate(smiles_list) if s not in completed_smiles]
    
    if not pending_smiles:
        print("\n✓ All SMILES already processed! Use --reset to rerun.")
        return []
    
    print(f"Pending: {len(pending_smiles)} (skipping {len(smiles_list) - len(pending_smiles)} already done)")
    
    print(f"\nSettings:")
    print(f"  LLM: {llm}")
    print(f"  AZ Model: {az_model}")
    print(f"  Hallucination Check: {hallucination_check} ({hallucination_method})")
    print(f"  Stability Check: {stability_check}")
    print(f"  Protecting Groups: {use_protecting_groups}")
    print(f"\n{'='*60}\n")
    
    results = []
    session_failed = []
    
    for progress_idx, (original_idx, smiles) in enumerate(pending_smiles):
        overall_done = len(completed_smiles) + progress_idx + 1
        print(f"[{overall_done}/{total_smiles}] Processing: {smiles[:50]}...")
        start_time = time.time()
        
        try:
            result = run_retrosynthesis(
                smiles=smiles,
                llm=llm,
                az_model=az_model,
                stability_flag=str(stability_check),
                hallucination_check=str(hallucination_check),
                use_protecting_group_feature=use_protecting_groups,
                hallucination_method=hallucination_method
            )
            
            elapsed = time.time() - start_time
            num_pathways = len(result.get('pathways', [])) if result else 0
            print(f"    ✓ Completed in {elapsed:.1f}s - {num_pathways} pathways found")
            
            # Save JSON result (same format as UI/API)
            json_path = save_result_json(smiles, result, results_subdir)
            print(f"    ✓ JSON saved to: {json_path.name}")
            
            # Update progress
            progress["completed"].append(smiles)
            progress["last_index"] = original_idx
            save_progress(input_file, progress)
            
            results.append({
                'smiles': smiles,
                'status': 'success',
                'num_pathways': num_pathways,
                'time': elapsed,
                'json_file': str(json_path),
                'result': result
            })
            
        except KeyboardInterrupt:
            print(f"\n\n*** INTERRUPTED - Progress saved! ***")
            print(f"    Completed: {len(progress['completed'])}/{total_smiles}")
            print(f"    Run again to resume from where you left off.")
            save_progress(input_file, progress)
            sys.exit(0)
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"    ✗ Failed in {elapsed:.1f}s - {str(e)[:50]}")
            
            # Still mark as "attempted" so we don't retry failures
            progress["failed"].append({"smiles": smiles, "error": str(e)})
            progress["completed"].append(smiles)  # Skip on retry
            save_progress(input_file, progress)
            
            session_failed.append({
                'smiles': smiles,
                'error': str(e)
            })
            results.append({
                'smiles': smiles,
                'status': 'failed',
                'error': str(e),
                'time': elapsed
            })
    
    # Summary
    total_completed = len(progress["completed"])
    total_failed = len(progress.get("failed", []))
    
    print(f"\n{'='*60}")
    print(f"BATCH COMPLETE")
    print(f"  Total in file: {total_smiles}")
    print(f"  Completed (all runs): {total_completed}")
    print(f"  Failed (all runs): {total_failed}")
    print(f"  This session: {len(pending_smiles)} processed")
    
    if session_failed:
        print(f"\nFailed this session:")
        for f in session_failed[:5]:
            print(f"  - {f['smiles'][:40]}... : {f['error'][:30]}")
        if len(session_failed) > 5:
            print(f"  ... and {len(session_failed)-5} more")
    
    print(f"\nLogs saved to: logs/{time.strftime('%Y-%m-%d')}/")
    print(f"Progress saved to: {get_progress_file(input_file)}")
    print(f"JSON results saved to: {results_subdir}")
    
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Batch retrosynthesis runner with auto-resume')
    parser.add_argument('--input', default='data/SMILES.txt', help='Input SMILES file')
    parser.add_argument('--llm', default='claude-sonnet-4-20250514', help='LLM model')
    parser.add_argument('--az-model', default='USPTO', help='AZ model')
    parser.add_argument('--hallucination-check', action='store_true', default=True, help='Enable hallucination check')
    parser.add_argument('--no-hallucination-check', action='store_false', dest='hallucination_check')
    parser.add_argument('--hallucination-method', default='rule_based', choices=['rule_based', 'ml_model'])
    parser.add_argument('--stability-check', action='store_true', help='Enable stability check')
    parser.add_argument('--protecting-groups', action='store_true', help='Enable protecting groups')
    parser.add_argument('--limit', type=int, help='Limit number of SMILES to process')
    parser.add_argument('--skip', type=int, default=0, help='Skip first N SMILES')
    parser.add_argument('--reset', action='store_true', help='Reset progress and start fresh')
    parser.add_argument('--clear-cache', action='store_true', help='Clear cache for all molecules before running')
    
    args = parser.parse_args()
    
    run_batch(
        input_file=args.input,
        llm=args.llm,
        az_model=args.az_model,
        stability_check=args.stability_check,
        hallucination_check=args.hallucination_check,
        hallucination_method=args.hallucination_method,
        use_protecting_groups=args.protecting_groups,
        limit=args.limit,
        skip=args.skip,
        reset=args.reset,
        clear_cache=args.clear_cache
    )

