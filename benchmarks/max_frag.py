"""
This script is used to compute the max fragment accuracy of the predictions.
It is based on the script from the following repository:
https://github.com/vthost/retroeval/blob/main/retroeval/utils/metrics.py
"""
from rdkit import Chem


def maxfrag(gt_reactants, pred_reactants):
    gt_reactants_mols = [
        Chem.MolFromSmiles(react) for react in gt_reactants.split(".")
    ]
    pred_reactants_mols = [
        Chem.MolFromSmiles(react) for react in pred_reactants.split(".")
    ]

    gt_reactants_maxfrag = Chem.MolToSmiles(
        max(gt_reactants_mols, key=lambda x: x.GetNumAtoms()))
    pred_reactants_maxfrag = Chem.MolToSmiles(
        max(pred_reactants_mols, key=lambda x: x.GetNumAtoms()))

    return gt_reactants_maxfrag, pred_reactants_maxfrag


def maxfrag_accuracy(gt_reactants, pred_reactants):
    gt_reactants_maxfrag, pred_reactants_maxfrag = maxfrag(
        gt_reactants, pred_reactants)
    return gt_reactants_maxfrag == pred_reactants_maxfrag
