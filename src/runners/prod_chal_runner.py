"""Runs retrosynthesis predictions on a predefined set of complex molecules.

This script is designed to process a hardcoded dictionary of challenging molecules
(defined in `mols_hard`). For each molecule in this set, it performs the following:
1. Retrieves the SMILES string.
2. Calls the `clear_cache()` function from `src.cache` before starting the loop.
3. Invokes the `main()` function (imported from `src.main`) with the molecule's
   SMILES string to get a result dictionary.
4. Saves this result dictionary as a JSON file in the
   `results/mols_hard/` directory. The filename includes the molecule's common
   name and a suffix (e.g., `_claude_cot.json`).

The script sets up the project root using `rootutils` to ensure correct
paths for imports and result saving.

Note:
    - The `mols_hard` dictionary contains several commented-out entries.
    - There is commented-out code that suggests an alternative processing path
      using a `run_prithvi` function, which is not currently active in this script.

Dependencies:
    - rootutils: For project root setup.
    - src.main: Provides the `main()` function for processing molecules.
    - src.cache: Provides the `clear_cache()` function.
    - json: For saving results.

Assumed File Structure:
    - A directory `results/mols_hard/` for saving output JSON files.
"""
# load up USPTO-50k test dataset
import rootutils

root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)
from src.main import main
from src.cache import clear_cache
import json

mols_hard = {
    # "maitotoxin":
    # "C[C@H](CC[C@@H]([C@@H]([C@H](C)C[C@H](C(=C)/C(=C/CO)/C)O)O)OS(=O)(=O)[O-])[C@H]([C@@H](C)[C@H]1[C@@H]([C@@H]([C@H]2[C@H](O1)[C@@H](C[C@]3([C@H](O2)C[C@H]4[C@H](O3)C[C@]5([C@H](O4)[C@H]([C@H]6[C@H](O5)C[C@H]([C@H](O6)[C@@H]([C@H](C[C@H]7[C@@H]([C@@H]([C@H]8[C@H](O7)C[C@H]9[C@H](O8)C[C@H]1[C@H](O9)[C@H]([C@@H]2[C@@H](O1)[C@@H]([C@H]([C@@H](O2)[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@H]1[C@@H](C[C@]2([C@H](O1)[C@@H]([C@]1([C@H](O2)C[C@]2([C@H](O1)CC[C@]1([C@H](O2)C[C@]2([C@H](O1)C[C@H]1[C@H](O2)CC[C@H](O1)[C@]1([C@@H](C[C@H]2[C@](O1)(C[C@H]1[C@](O2)(CC[C@]2([C@H](O1)C[C@H]1[C@](O2)(C[C@H]2[C@H](O1)C/C=C\[C@H]1[C@H](O2)C[C@H]2[C@](O1)(C[C@]1([C@H](O2)C[C@H]2[C@](O1)(CC[C@H](O2)[C@H]([C@@H](C[C@@H](C)[C@@H](C)CC=C)O)O)C)C)C)C)C)C)C)O)C)C)C)C)C)O)C)O)O)O)O)O)O)O)O)O)O)O)O)O)OS(=O)(=O)[O-])O)O)O)O)C)C)O)O)O)O",
    # "taxol":
    # "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@](C3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
    "dalbavancin":
    "CC(C)CCCCCCCCC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1OC2=C3C=C4C=C2OC5=C(C=C(C=C5)[C@H]([C@H]6C(=O)N[C@@H](C7=C(C(=CC(=C7)O)O[C@@H]8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)C9=C(C=CC(=C9)[C@H](C(=O)N6)NC(=O)[C@@H]4NC(=O)[C@@H]1C2=C(C(=CC(=C2)OC2=C(C=CC(=C2)[C@H](C(=O)N[C@H](CC2=CC=C(O3)C=C2)C(=O)N1)NC)O)O)Cl)O)C(=O)NCCCN(C)C)O)Cl)C(=O)O)O)O",
    # "palytoxin":
    # "C[C@H]1C[C@@]2([C@H](O[C@](C1)(O2)CCCCCCC[C@@H](C[C@@H]3[C@@H]([C@H]([C@H]([C@@](O3)(C[C@@H]([C@@H](C)/C=C/[C@H](CC[C@H]([C@H]([C@@H]4C[C@H]([C@@H]([C@H](O4)C[C@H]([C@@H](C[C@@H]5[C@H]([C@@H]([C@H]([C@@H](O5)C[C@@H](/C=C\C=C\C[C@H]([C@@H]([C@@H](C/C=C\C(=C)CC[C@@H]([C@H]([C@@H]([C@H](C)C[C@@H]6[C@@H]([C@H]([C@@H]([C@H](O6)/C=C\[C@H]([C@@H](C[C@@H]7C[C@@H]8C[C@H](O7)[C@H](O8)CC[C@@H]9[C@@H](C[C@H](O9)CN)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)C[C@@H](C)CCCCC[C@H]([C@@H]([C@@H]([C@H]([C@@H]([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)C[C@@H]([C@@H](/C(=C/[C@@H](C[C@@H](C)[C@@H](C(=O)N/C=C/C(=O)NCCCO)O)O)/C)O)O)O)O)O)O)O)O)O)O)C",
    # "Sirolimus":
    # "C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC",
    "cyclosporin a":
    "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C",
    # "rapamycin":
    # "C[C@@H]1CCC2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC",
}

# # now run prithvi on the dataset
clear_cache()

for mol in mols_hard:
    molecule = mols_hard[mol]
    print(f"Running {mol}")
    # try:
    #     res_dict = run_prithvi(molecule)
    #     with open(f"{root_dir}/results/mols_hard/{mol}.json", "w") as f:
    #         json.dump(res_dict, f, indent=4)
    # except Exception as e:
    #     print("Error in molecule:", mol)
    #     print("Error:", e)
    res_dict = main(molecule)  #, llm="o1-preview-2024-09-12")
    with open(f"{root_dir}/results/mols_hard/{mol}_claude_cot.json", "w") as f:
        json.dump(res_dict, f, indent=4)
