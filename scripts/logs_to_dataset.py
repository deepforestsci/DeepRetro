"""
Convert log files to ML training dataset.

Usage:
    python scripts/logs_to_dataset.py                           # Process all logs
    python scripts/logs_to_dataset.py --date 2025-12-16         # Specific date
    python scripts/logs_to_dataset.py --job job_20251216_084014 # Specific job
    python scripts/logs_to_dataset.py --output csv              # Output as CSV
    python scripts/logs_to_dataset.py --output jsonl            # Output as JSONL
"""

import os
import re
import ast
import json
import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime


def parse_hallucination_log_line(line: str) -> dict:
    """Parse a hallucination check log line into structured data."""
    
    # Pattern to match: Hallucination check - pathway_idx: X, product: X, reactants: X, filtered_out: X, report: {...}
    pattern = r"Hallucination check - pathway_idx: (\d+), product: ([^,]+), reactants: ([^,]+), filtered_out: (True|False), report: (\{.*\})"
    
    match = re.search(pattern, line)
    if not match:
        return None
    
    try:
        report = ast.literal_eval(match.group(5))
    except:
        return None
    
    # Extract molecular features if present
    details = report.get('details', {})
    mol_features = details.get('molecular_features', {})
    
    record = {
        # Timestamp (extract from log line start)
        'timestamp': extract_timestamp(line),
        
        # Core fields
        'pathway_idx': int(match.group(1)),
        'product_smiles': match.group(2).strip(),
        'reactants_smiles': match.group(3).strip(),
        
        # Label (for ML training)
        'filtered_out': match.group(4) == 'True',
        'label': 1 if match.group(4) == 'True' else 0,
        
        # Score
        'score': report.get('score'),
        'severity': report.get('severity'),
        
        # Comparison results
        'valid_reactant': details.get('valid_reactant'),
        'valid_product': details.get('valid_product'),
        'atom_count_consistent': details.get('atom_count_consistent'),
        'num_ring_changes': len(details.get('ring_size_changes', [])),
        'num_substituent_changes': len(details.get('substituent_position_changes', [])),
        'num_issues': len(details.get('detected_issues', [])),
        
        # Molecular features (flattened for ML)
        'reactant_num_aromatic': mol_features.get('reactant_num_aromatic'),
        'product_num_aromatic': mol_features.get('product_num_aromatic'),
        'reactant_total_bonds': mol_features.get('reactant_total_bonds'),
        'product_total_bonds': mol_features.get('product_total_bonds'),
        'reactant_ring_count': len(mol_features.get('reactant_ring_sizes', [])),
        'product_ring_count': len(mol_features.get('product_ring_sizes', [])),
        
        # Raw data (for detailed analysis)
        'penalties': report.get('penalties', []),
        'detected_issues': details.get('detected_issues', []),
        'ring_size_changes': details.get('ring_size_changes', []),
        'reactant_atom_counts': mol_features.get('reactant_atom_counts', {}),
        'product_atom_counts': mol_features.get('product_atom_counts', {}),
    }
    
    # Calculate derived features
    if record['reactant_atom_counts'] and record['product_atom_counts']:
        record['atom_count_diff'] = sum(record['reactant_atom_counts'].values()) - sum(record['product_atom_counts'].values())
    else:
        record['atom_count_diff'] = None
        
    if record['reactant_total_bonds'] and record['product_total_bonds']:
        record['bond_count_diff'] = record['reactant_total_bonds'] - record['product_total_bonds']
    else:
        record['bond_count_diff'] = None
    
    return record


def parse_stability_log_line(line: str) -> dict:
    """Parse a stability check log line into structured data."""
    
    # Pattern: Stability dict: {...}
    pattern = r"Stability dict: (\{.*\})"
    
    match = re.search(pattern, line)
    if not match:
        return None
    
    try:
        stability_dict = ast.literal_eval(match.group(1))
    except:
        return None
    
    metrics = stability_dict.get('metrics', {})
    ring_data = stability_dict.get('ring_data', {})
    atom_data = stability_dict.get('atom_data', {})
    
    # Determine if filtered (assessment not in stable categories)
    assessment = stability_dict.get('assessment', '')
    filtered_out = assessment not in ['Likely stable', 'Moderately stable']
    
    record = {
        'timestamp': extract_timestamp(line),
        'check_type': 'stability',
        
        # Label
        'filtered_out': filtered_out,
        'label': 1 if filtered_out else 0,
        
        # Score
        'stability_score': stability_dict.get('stability_score'),
        'assessment': assessment,
        'num_issues': len(stability_dict.get('issues', [])),
        
        # Metrics
        'molecular_weight': metrics.get('molecular_weight'),
        'logP': metrics.get('logP'),
        'h_bond_donors': metrics.get('h_bond_donors'),
        'h_bond_acceptors': metrics.get('h_bond_acceptors'),
        'rotatable_bonds': metrics.get('rotatable_bonds'),
        
        # Ring data
        'num_rings': ring_data.get('num_rings'),
        'num_aliphatic_rings': ring_data.get('num_aliphatic_rings'),
        'num_aromatic_rings': ring_data.get('num_aromatic_rings'),
        'num_bridgehead_atoms': ring_data.get('num_bridgehead_atoms'),
        
        # Atom data
        'num_atoms': atom_data.get('num_atoms'),
        'num_bonds': atom_data.get('num_bonds'),
        'num_heavy_atoms': atom_data.get('num_heavy_atoms'),
        'num_aromatic_atoms': atom_data.get('num_aromatic_atoms'),
        
        # Raw
        'issues': stability_dict.get('issues', []),
    }
    
    return record


def extract_timestamp(line: str) -> str:
    """Extract timestamp from log line."""
    # Pattern: 2025-12-16 03:10:14
    match = re.match(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)
    if match:
        return match.group(1)
    return None


def extract_job_id(line: str) -> str:
    """Extract job_id from log line."""
    match = re.search(r'job_id=(\S+)', line)
    if match:
        return match.group(1)
    return None


def process_log_file(filepath: str) -> tuple:
    """Process a single log file and extract all check records."""
    
    hallucination_records = []
    stability_records = []
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # Extract job_id for context
            job_id = extract_job_id(line)
            
            # Try hallucination check
            if 'Hallucination check -' in line:
                record = parse_hallucination_log_line(line)
                if record:
                    record['job_id'] = job_id
                    record['source_file'] = os.path.basename(filepath)
                    hallucination_records.append(record)
            
            # Try stability check
            elif 'Stability dict:' in line:
                record = parse_stability_log_line(line)
                if record:
                    record['job_id'] = job_id
                    record['source_file'] = os.path.basename(filepath)
                    stability_records.append(record)
    
    return hallucination_records, stability_records


def find_log_files(logs_dir: str, date: str = None, job: str = None) -> list:
    """Find log files matching criteria."""
    
    log_files = []
    logs_path = Path(logs_dir)
    
    if not logs_path.exists():
        print(f"Logs directory not found: {logs_dir}")
        return []
    
    if job:
        # Find specific job file
        for f in logs_path.rglob(f"*{job}*.log"):
            log_files.append(str(f))
    elif date:
        # Find all files for a specific date
        date_dir = logs_path / date
        if date_dir.exists():
            for f in date_dir.glob("*.log"):
                log_files.append(str(f))
    else:
        # Find all log files
        for f in logs_path.rglob("*.log"):
            log_files.append(str(f))
    
    return sorted(log_files)


def save_dataset(records: list, output_path: str, format: str = 'csv'):
    """Save records to file."""
    
    if not records:
        print(f"No records to save")
        return
    
    if format == 'csv':
        # For CSV, we need to flatten list/dict columns
        df = pd.DataFrame(records)
        
        # Convert list columns to string for CSV
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].apply(lambda x: json.dumps(x) if isinstance(x, (list, dict)) else x)
        
        df.to_csv(output_path, index=False)
        print(f"Saved {len(records)} records to {output_path}")
        
    elif format == 'jsonl':
        with open(output_path, 'w') as f:
            for record in records:
                f.write(json.dumps(record) + '\n')
        print(f"Saved {len(records)} records to {output_path}")
        
    elif format == 'json':
        with open(output_path, 'w') as f:
            json.dump(records, f, indent=2)
        print(f"Saved {len(records)} records to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Convert logs to ML training dataset')
    parser.add_argument('--logs-dir', default='logs', help='Logs directory path')
    parser.add_argument('--date', help='Process logs for specific date (YYYY-MM-DD)')
    parser.add_argument('--job', help='Process logs for specific job ID')
    parser.add_argument('--output', default='csv', choices=['csv', 'jsonl', 'json'], help='Output format')
    parser.add_argument('--output-dir', default='training_data', help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Find log files
    log_files = find_log_files(args.logs_dir, args.date, args.job)
    print(f"Found {len(log_files)} log files")
    
    if not log_files:
        return
    
    # Process all files
    all_hallucination = []
    all_stability = []
    
    for filepath in log_files:
        print(f"Processing: {filepath}")
        h_records, s_records = process_log_file(filepath)
        all_hallucination.extend(h_records)
        all_stability.extend(s_records)
    
    # Generate timestamp for output files
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Save datasets
    if all_hallucination:
        output_file = f"{args.output_dir}/hallucination_dataset_{timestamp}.{args.output}"
        save_dataset(all_hallucination, output_file, args.output)
        
        # Print summary
        df = pd.DataFrame(all_hallucination)
        print(f"\nHallucination Dataset Summary:")
        print(f"  Total records: {len(df)}")
        print(f"  Filtered out (label=1): {df['label'].sum()}")
        print(f"  Kept (label=0): {len(df) - df['label'].sum()}")
        if 'severity' in df.columns:
            print(f"  Severity distribution: {df['severity'].value_counts().to_dict()}")
    
    if all_stability:
        output_file = f"{args.output_dir}/stability_dataset_{timestamp}.{args.output}"
        save_dataset(all_stability, output_file, args.output)
        
        # Print summary
        df = pd.DataFrame(all_stability)
        print(f"\nStability Dataset Summary:")
        print(f"  Total records: {len(df)}")
        print(f"  Filtered out (label=1): {df['label'].sum()}")
        print(f"  Kept (label=0): {len(df) - df['label'].sum()}")
        if 'assessment' in df.columns:
            print(f"  Assessment distribution: {df['assessment'].value_counts().to_dict()}")


if __name__ == '__main__':
    main()

