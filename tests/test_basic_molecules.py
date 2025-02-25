import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

from src.utils.llm import call_LLM, split_cot_json, split_json_openAI, split_json_deepseek
from src.variables import OPENAI_MODELS, DEEPSEEK_MODELS

claude_adv = "claude-3-opus-20240229:adv"

openai_adv = "gpt-4o:adv"

deepseek_adv = "azure_ai/DeepSeek-R1:adv"

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

def test_basic_molecules_adv_prompt_claude(llm=claude_adv):
    from src.variables import BASIC_MOLECULES

    print("BASIC_MOLECULES: ", BASIC_MOLECULES)
    successful_molecules = []
    failed_molecules = []

    for molecule in BASIC_MOLECULES:
        try:
            status_code, res = call_LLM(molecule=molecule, LLM = llm)
            assert status_code == 200
            
            if not res:
                raise Exception("Response is empty.")
            
            successful_molecules.append(molecule)
        except Exception as e:
            print(f"Error: {e}\n Molecule: {molecule}")
            failed_molecules.append(molecule)
    if failed_molecules:
        print(f"Failed molecules: {failed_molecules}")

def test_basic_molecules_adv_prompt_openai(llm=openai_adv): 
    from src.variables import BASIC_MOLECULES

    print("BASIC_MOLECULES: ", BASIC_MOLECULES)
    successful_molecules = []
    failed_molecules = []

    for molecule in BASIC_MOLECULES:
        try:
            status_code, res = call_LLM(molecule=molecule, LLM = llm)
            assert status_code == 200
            
            if not res:
                raise Exception("Response is empty.")
            
            successful_molecules.append(molecule)
        except Exception as e:
            print(f"Error: {e}\n Molecule: {molecule}")
            failed_molecules.append(molecule)
    if failed_molecules:
        print(f"Failed molecules: {failed_molecules}")

def test_basic_molecules_adv_prompt_deepseek(llm=deepseek_adv):
    from src.variables import BASIC_MOLECULES

    print("BASIC_MOLECULES: ", BASIC_MOLECULES)
    successful_molecules = []
    failed_molecules = []

    for molecule in BASIC_MOLECULES:
        try:
            status_code, _ = call_LLM(molecule=molecule, LLM = llm)
            assert status_code == 200
            successful_molecules.append(molecule)
        except Exception as e:
            print(f"Error: {e}\n Molecule: {molecule}")
            failed_molecules.append(molecule)
    if failed_molecules:
        print(f"Failed molecules: {failed_molecules}")

if __name__ == '__main__':
    pytest.main()