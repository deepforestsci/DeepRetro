import unittest
from unittest.mock import patch, ANY
from src.main import main
from src.utils.job_context import logger as context_logger

class TestMain(unittest.TestCase):
    def setUp(self):
        self.test_smiles = "CCO"
        self.test_llm = "claude-3-opus-20240229"
        self.test_az_model = "USPTO"
        self.test_stability_flag = "False"
        self.test_hallucination_check = "False"

    @patch('src.main.run_prithvi')
    @patch('src.main.setup_logging')
    def test_main_success(self, mock_setup_logging, mock_run_prithvi):
        # Mock successful run
        mock_run_prithvi.return_value = {"result": "test_result"}
        
        result = main(
            smiles=self.test_smiles,
            llm=self.test_llm,
            az_model=self.test_az_model,
            stability_flag=self.test_stability_flag,
            hallucination_check=self.test_hallucination_check
        )
        
        self.assertEqual(result, {"result": "test_result"})
        mock_setup_logging.assert_called_once()
        mock_run_prithvi.assert_called_once_with(
            molecule=self.test_smiles,
            llm=self.test_llm,
            az_model=self.test_az_model,
            stability_flag=self.test_stability_flag,
            hallucination_check=self.test_hallucination_check
        )

    @patch('src.main.run_prithvi')
    @patch('src.main.setup_logging')
    def test_main_run_prithvi_exception(self, mock_setup_logging, mock_run_prithvi):
        # Mock run_prithvi raising an exception
        mock_run_prithvi.side_effect = Exception("Test error")
        
        with self.assertRaises(Exception) as context:
            main(
                smiles=self.test_smiles,
                llm=self.test_llm,
                az_model=self.test_az_model,
                stability_flag=self.test_stability_flag,
                hallucination_check=self.test_hallucination_check
            )
        
        self.assertEqual(str(context.exception), "Test error")
        mock_setup_logging.assert_called_once()
        mock_run_prithvi.assert_called_once_with(
            molecule=self.test_smiles,
            llm=self.test_llm,
            az_model=self.test_az_model,
            stability_flag=self.test_stability_flag,
            hallucination_check=self.test_hallucination_check
        )

    @patch('src.main.run_prithvi')
    @patch('src.main.setup_logging')
    def test_main_setup_logging_exception(self, mock_setup_logging, mock_run_prithvi):
        # Mock setup_logging raising an exception
        mock_setup_logging.side_effect = Exception("Logging setup error")
        
        with self.assertRaises(Exception) as context:
            main(
                smiles=self.test_smiles,
                llm=self.test_llm,
                az_model=self.test_az_model,
                stability_flag=self.test_stability_flag,
                hallucination_check=self.test_hallucination_check
            )
        
        self.assertEqual(str(context.exception), "Logging setup error")
        mock_setup_logging.assert_called_once()
        mock_run_prithvi.assert_not_called()

if __name__ == '__main__':
    unittest.main() 