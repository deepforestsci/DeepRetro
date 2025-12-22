import json
import csv
from pathlib import Path

INPUT_PATH = Path("data/checker_dataset_pistachio.jsonl")
OUTPUT_PATH = Path("data/checker_dataset_clean.csv")

def main():
    rows = []

    with INPUT_PATH.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)

            target = rec.get("target_smiles_canonical")
            product = rec.get("product_smiles_canonical")
            hall = rec.get("hallucination_raw", {}) or {}

            hallucination_score = hall.get("score")
            hallucination_severity = hall.get("severity")

            # stability_raw is a dict: reactant_smiles -> stability_dict
            stability_raw = rec.get("stability_raw", {}) or {}

            # reactant_smiles_canonical is a list aligned with reactant_smiles
            reactants_canonical = rec.get("reactant_smiles_canonical", [])

            for rcanon in reactants_canonical:
                stab = stability_raw.get(rcanon, {})  # prefer canonical key if present
                if not stab:
                    # fallback: try by non-canonical reactant name if needed
                    for key, val in stability_raw.items():
                        if key == rcanon:
                            stab = val
                            break

                stability_score = stab.get("stability_score")
                stability_assessment = stab.get("assessment")
                issues = stab.get("issues", []) or []
                stability_issues = "; ".join(issues)

                rows.append({
                    "target_smiles_canonical": target,
                    "product_smiles_canonical": product,
                    "reactant_smiles_canonical": rcanon,
                    "stability_score": stability_score,
                    "stability_assessment": stability_assessment,
                    "stability_issues": stability_issues,
                    "hallucination_score": hallucination_score,
                    "hallucination_severity": hallucination_severity,
                })

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", encoding="utf-8", newline="") as out_f:
        writer = csv.DictWriter(
            out_f,
            fieldnames=[
                "target_smiles_canonical",
                "product_smiles_canonical",
                "reactant_smiles_canonical",
                "stability_score",
                "stability_assessment",
                "stability_issues",
                "hallucination_score",
                "hallucination_severity",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()