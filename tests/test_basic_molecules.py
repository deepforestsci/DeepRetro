import os
import ast
import pytest

import rootutils
root_dir = rootutils.setup_root(".", indicator=".project-root", pythonpath=True)

from src.utils.llm import call_LLM

def test_basic_molecules():
    from src.variables import BASIC_MOLECULES

    for molecule in BASIC_MOLECULES:
        try:
            status_code, _ = call_LLM(molecule=molecule)
            assert status_code == 200
        except Exception as e:
            print(f"Error: {e}\n Molecule: {molecule}")

if __name__ == '__main__':
    pytest.main()