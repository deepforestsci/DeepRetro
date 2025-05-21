import unittest
from unittest.mock import patch, ANY
from src.metadata import reagent_agent, conditions_agent, literature_agent
from src.utils.job_context import logger as context_logger

class TestMetadata(unittest.TestCase):
    def setUp(self):
        self.test_reactants = [{"smiles": "CCO", "name": "ethanol"}]
        self.test_products = [{"smiles": "CC(=O)O", "name": "acetic acid"}]
        self.test_reagents = [{"smiles": "O", "name": "water"}]
        self.test_conditions = {
            "temperature": "25°C",
            "pressure": "1 atm",
            "solvent": "water",
            "time": "1 hour"
        }

if __name__ == '__main__':
    unittest.main() 