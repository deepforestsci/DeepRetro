import rootutils

root_dir = rootutils.setup_root(".",
                                indicator=".project-root",
                                pythonpath=True)

import pandas as pd
import json
from rdkit import Chem


def get_results_in_file(file):
    step_1 = file['steps'][0]

    reactants = []
    for react in step_1['reactants']:
        reactants.append(react['smiles'])
    reagents = []
    for reagent in step_1['reagents']:
        reagents.append(reagent['smiles'])
    products = []
    for prod in step_1['products']:
        products.append(prod['smiles'])
    return reactants, reagents, products


def all_correct(input, output) -> bool:
    input_splits = input.split('.')
    output_splits = output.split('.')
    input_mols = [Chem.CanonSmiles(mol) for mol in input_splits]
    output_mols = [Chem.CanonSmiles(mol) for mol in output_splits]
    if all([output_mol in input_mols for output_mol in output_mols]):
        return 1
    else:
        return 0


def any_correct(input, output) -> bool:
    input_splits = input.split('.')
    output_splits = output.split('.')
    input_mols = [Chem.CanonSmiles(mol) for mol in input_splits]
    output_mols = [Chem.CanonSmiles(mol) for mol in output_splits]
    if any([output_mol in input_mols for output_mol in output_mols]):
        return 1
    else:
        return 0


def new_eval(df_results):
    df_results['all_correct'] = df_results.apply(
        lambda row: all_correct(row['input'], row['output']), axis=1)
    df_results['any_correct'] = df_results.apply(
        lambda row: any_correct(row['input'], row['output']), axis=1)
    all_correct = 100 * df_results.all_correct.sum() / len(df_results)
    any_correct = 100 * df_results.any_correct.sum() / len(df_results)
    return all_correct, any_correct, df_results


def get_dataset(path, canonicalize=True):
    df = pd.read_csv(path)
    # canonicalize the input
    if canonicalize:
        df['input'] = df['input'].apply(Chem.CanonSmiles)
    return df


def get_results_csv(path):
    df = pd.read_csv(path)
    return df


def merge_results(df, results_df):
    # merge the results_df with the df based on the input column
    df = df.join(results_df, on='input')
    return df


def main():
    df = get_dataset("../benchmarks/sampled_molecules.csv")


def actual_eval(df_results):
    # check if the predicted outputs are in the actual outputs
    # and if the number of predicted outputs is equal to the number of actual outputs
    not_common = []
    missed = []
    for i, row in df_results.iterrows():
        outputs = row['output'].split('.')
        # print(outputs)
        predicted_outputs = row['reactants']  # + row['reagents']
        # print(outputs, predicted_outputs)
        outputs_mols = [Chem.CanonSmiles(output) for output in outputs]
        predicted_outputs_mols = [
            Chem.CanonSmiles(output) for output in predicted_outputs
        ]
        # print(outputs_mols,predicted_outputs_mols)
        # break
        # check if all the predicted outputs are in the actual outputs
        if any([output in outputs_mols for output in predicted_outputs_mols]):
            df_results.loc[i, 'any_correct'] = 1
        else:
            df_results.loc[i, 'any_correct'] = 0

        if all([output in outputs_mols for output in predicted_outputs_mols]):
            df_results.loc[i, 'all_correct'] = 1
        else:
            df_results.loc[i, 'all_correct'] = 0
        # get the molecules that are not common
        not_common.append([
            output for output in predicted_outputs_mols
            if output not in outputs_mols
        ])
        # df_results.loc[i, 'not_common'] =
        # get the molecules that are missed
        missed.append([
            output for output in outputs_mols
            if output not in predicted_outputs_mols
        ])

    df_results['not_common'] = not_common
    df_results['missed'] = missed
    all_correct = 100 * df_results.all_correct.sum() / len(df_results)
    any_correct = 100 * df_results.any_correct.sum() / len(df_results)
    return all_correct, any_correct, df_results
