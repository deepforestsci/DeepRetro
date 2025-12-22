"""Offline runner to collect rule-based checker data over a SMILES dataset.

This script runs multi-step retrosynthesis and collects outputs from
rule-based validation engines (stability and hallucination checkers)
for training purposes. It does not modify any core pipeline code.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any, Dict, List, Tuple, Optional

from rdkit import Chem
import structlog

from src.rec_prithvi import rec_run_prithvi
from src.utils.stability_checks import check_molecule_stability
from src.utils.hallucination_checks import calculate_hallucination_score
from src.utils.job_context import logger as context_logger


def _canonical_smiles(smiles: str) -> str:
    """Return canonical SMILES string using RDKit.

    Parameters
    ----------
    smiles : str
        Input SMILES string

    Returns
    -------
    str
        Canonical SMILES, or original string if parsing fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return smiles


def _read_smiles_file(path: str) -> List[Tuple[str, Optional[str]]]:
    """Read SMILES dataset file.

    Parameters
    ----------
    path : str
        Path to input file (one SMILES per line, optional ID)

    Returns
    -------
    List[Tuple[str, Optional[str]]]
        List of (smiles, optional_id) tuples
    """
    entries: List[Tuple[str, Optional[str]]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            mol_id = parts[1] if len(parts) > 1 else None
            entries.append((smiles, mol_id))
    return entries


def _traverse_tree_and_collect_steps(
    root: Dict[str, Any],
) -> List[Tuple[str, List[str], int, Optional[int], int]]:
    """Traverse retrosynthesis tree and extract reaction steps.

    Parameters
    ----------
    root : Dict[str, Any]
        Root node from rec_run_prithvi result

    Returns
    -------
    List[Tuple[str, List[str], int, Optional[int], int]]
        List of (product_smiles, reactant_smiles_list, step_id, parent_step_id, depth)
    """
    steps: List[Tuple[str, List[str], int, Optional[int], int]] = []

    step_counter = 0

    def visit_mol_node(
        mol_node: Dict[str, Any],
        parent_step_id: Optional[int],
        depth: int,
    ) -> None:
        nonlocal step_counter
        product_smiles = mol_node.get("smiles", "")
        children = mol_node.get("children") or []

        for child in children:
            if child.get("type") != "reaction":
                continue

            reactant_smiles: List[str] = []
            for reactant_node in child.get("children") or []:
                if reactant_node.get("type") == "mol":
                    reactant_smiles.append(reactant_node.get("smiles", ""))

            step_id = step_counter
            step_counter += 1

            steps.append(
                (product_smiles, reactant_smiles, step_id, parent_step_id, depth)
            )

            for reactant_node in child.get("children") or []:
                if reactant_node.get("type") == "mol":
                    visit_mol_node(
                        reactant_node,
                        parent_step_id=step_id,
                        depth=depth + 1,
                    )

    if root.get("type") == "mol":
        visit_mol_node(root, parent_step_id=None, depth=0)

    return steps


def collect_checker_data_for_target(
    target_smiles: str,
    target_id: Optional[str],
    llm: str,
    az_model: str,
) -> List[Dict[str, Any]]:
    """Run retrosynthesis and collect checker outputs for a target molecule.

    Parameters
    ----------
    target_smiles : str
        Target molecule SMILES
    target_id : Optional[str]
        Optional identifier for the target
    llm : str
        LLM model name
    az_model : str
        AiZynthFinder model name

    Returns
    -------
    List[Dict[str, Any]]
        List of records, one per reaction step
    """
    job_id = f"checker_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"

    log = structlog.get_logger().bind(job_id=job_id)
    token = context_logger.set(log)

    try:
        print(f"[checker] Starting job_id={job_id} for target={target_smiles!r}")

        # Run the core retrosynthesis (multi-step)
        result_dict, solved = rec_run_prithvi(
            molecule=target_smiles,
            job_id=job_id,
            llm=llm,
            az_model=az_model,
            stability_flag="True",
            hallucination_check="True",
            use_protecting_group_feature=False,
        )
    finally:
        context_logger.reset(token)

    steps = _traverse_tree_and_collect_steps(result_dict)

    records: List[Dict[str, Any]] = []
    target_canonical = _canonical_smiles(target_smiles)
    timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")

    print(
        f"[checker] Finished retrosynthesis for target={target_smiles!r}, "
        f"solved={bool(solved)}, steps_found={len(steps)}"
    )

    for product_smiles, reactant_smiles, step_id, parent_step_id, depth in steps:
        stability_results: Dict[str, Any] = {}
        for rs in reactant_smiles:
            stability_results[rs] = check_molecule_stability(rs)

        combined_reactants = ".".join(reactant_smiles)
        hallucination_result = calculate_hallucination_score(
            reactant_smiles=combined_reactants,
            product_smiles=product_smiles,
        )

        record: Dict[str, Any] = {
            "timestamp": timestamp,
            "job_id": job_id,
            "target_smiles_raw": target_smiles,
            "target_smiles_canonical": target_canonical,
            "target_id": target_id,
            "llm": llm,
            "az_model": az_model,
            "solved": bool(solved),
            "step_id": step_id,
            "parent_step_id": parent_step_id,
            "depth": depth,
            "product_smiles": product_smiles,
            "product_smiles_canonical": _canonical_smiles(product_smiles),
            "reactant_smiles": reactant_smiles,
            "reactant_smiles_canonical": [
                _canonical_smiles(rs) for rs in reactant_smiles
            ],
            "stability_raw": stability_results,
            "hallucination_raw": hallucination_result,
        }
        records.append(record)

    return records


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect rule-based checker data from retrosynthesis runs"
    )
    parser.add_argument(
        "--input",
        default="data/targets.smi",
        help="Path to input SMILES file (one SMILES per line, optional ID)",
    )
    parser.add_argument(
        "--output",
        default="data/checker_dataset_pistachio.jsonl",
        help="Path to output JSONL file",
    )
    parser.add_argument(
        "--llm",
        default="claude-opus-4-20250514",
        help="LLM model name",
    )
    parser.add_argument(
        "--az_model",
        default="USPTO",
        help="AiZynthFinder model name",
    )

    args = parser.parse_args()

    entries = _read_smiles_file(args.input)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    import gc

    with open(args.output, "a", encoding="utf-8") as out_f:
        for idx, (smiles, mol_id) in enumerate(entries):
            records = collect_checker_data_for_target(
                target_smiles=smiles,
                target_id=mol_id,
                llm=args.llm,
                az_model=args.az_model,
            )
            for rec in records:
                out_f.write(json.dumps(rec) + "\n")
            out_f.flush()
            del records
            gc.collect()


if __name__ == "__main__":
    main()


