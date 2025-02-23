'''
Note: This set of tests takes quite some time to run, as it runs over 86 BASIC_MOLECULES listed in the src/variables.py
'''

import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

from src.utils.llm import call_LLM

def test_basic_molecules():
    from src.variables import BASIC_MOLECULES

    print("BASIC_MOLECULES: ", BASIC_MOLECULES)
    successful_molecules = []
    failed_molecules = []

    for molecule in BASIC_MOLECULES:
        try:
            status_code, res = call_LLM(molecule=molecule)
            assert status_code == 200

            if not res:
                raise Exception("Response is empty.")
            
            successful_molecules.append(molecule)

        except Exception as e:
            print(f"Error: {e}\n Molecule: {molecule}")
            failed_molecules.append(molecule)
    if failed_molecules:
        print(f"Failed molecules: {failed_molecules}")

# def test_basic_molecules_adv_prompt_openai(): 
#     from src.variables import BASIC_MOLECULES
#     from tests.variables_test import OPENAI_ADV_MODEL

#     print("BASIC_MOLECULES: ", BASIC_MOLECULES)
#     successful_molecules = []
#     failed_molecules = []

#     for molecule in BASIC_MOLECULES:
#         try:
#             status_code, res = call_LLM(molecule=molecule, LLM = OPENAI_ADV_MODEL)
#             assert status_code == 200
            
#             if not res:
#                 raise Exception("Response is empty.")
            
#             successful_molecules.append(molecule)
#         except Exception as e:
#             print(f"Error: {e}\n Molecule: {molecule}")
#             failed_molecules.append(molecule)
#     if failed_molecules:
#         print(f"Failed molecules: {failed_molecules}")

def test_basic_molecules_adv_prompt_deepseek():
    from src.variables import BASIC_MOLECULES
    from tests.variables_test import DEEPSEEK_FIREWORKS_MODEL

    print("BASIC_MOLECULES: ", BASIC_MOLECULES)
    successful_molecules = []
    failed_molecules = []

    for molecule in BASIC_MOLECULES:
        try:
            status_code, _ = call_LLM(molecule=molecule,
                                      LLM = DEEPSEEK_FIREWORKS_MODEL)
            assert status_code == 200
            successful_molecules.append(molecule)
        except Exception as e:
            print(f"Error: {e}\n Molecule: {molecule}")
            failed_molecules.append(molecule)
    if failed_molecules:
        print(f"Failed molecules: {failed_molecules}")

if __name__ == '__main__':
    pytest.main()