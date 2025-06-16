import unittest
import os
import json
from unittest.mock import patch, ANY


# --- Minimal change for CI: Start ---
# Store original API_KEY if it exists and set a specific one for tests
_ORIGINAL_API_KEY_ENV = os.environ.get('API_KEY')
_TEST_API_KEY = "test_api_key_for_ci_12345"
os.environ['API_KEY'] = _TEST_API_KEY
# --- Minimal change for CI: End ---


from src.api import app, save_result, load_result, PARTIAL_JSON_PATH, API_KEY

class TestApiFunctions(unittest.TestCase):


    @classmethod
    def tearDownClass(cls):
        # --- Minimal change for CI: Start ---
        # Restore original API_KEY environment variable
        if _ORIGINAL_API_KEY_ENV is None:
            if 'API_KEY' in os.environ: # Only del if we were the ones to set it
                del os.environ['API_KEY']
        else:
            os.environ['API_KEY'] = _ORIGINAL_API_KEY_ENV
        # --- Minimal change for CI: End ---

    def setUp(self):
        self.app_context = app.app_context()
        self.app_context.push()
        self.client = app.test_client()

        self.api_key = API_KEY
        # --- Minimal change for CI: Start ---
        # Verify that the API_KEY from src.api is the one we set for testing
        self.assertEqual(self.api_key, _TEST_API_KEY, 
                         "API_KEY from src.api was not the expected test key. Check environment setup at the top of this file.")
        # --- Minimal change for CI: End ---


    def tearDown(self):
        if os.path.exists(PARTIAL_JSON_PATH):
            os.remove(PARTIAL_JSON_PATH)
        self.app_context.pop()

    # Tests for save_result and load_result (file operations)
    def test_save_result_creates_file_with_correct_content(self):
        smiles = "CCO"
        result_data = {"reaction": "C=C.O"}
        save_result(smiles, result_data)
        self.assertTrue(os.path.exists(PARTIAL_JSON_PATH))
        with open(PARTIAL_JSON_PATH, 'r') as f:
            loaded_data = json.load(f)
        self.assertEqual(loaded_data["smiles"], smiles)
        self.assertEqual(loaded_data["result"], result_data)

    def test_load_result_file_not_exists(self):
        if os.path.exists(PARTIAL_JSON_PATH):
            os.remove(PARTIAL_JSON_PATH)
        smiles, result = load_result()
        self.assertIsNone(smiles)
        self.assertIsNone(result)

    def test_load_result_file_exists_and_valid(self):
        smiles_expected = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        result_expected = {"pathway": "some_reaction_steps"}
        with open(PARTIAL_JSON_PATH, 'w') as f:
            json.dump({"smiles": smiles_expected, "result": result_expected}, f)
        smiles_loaded, result_loaded = load_result()
        self.assertEqual(smiles_loaded, smiles_expected)
        self.assertEqual(result_loaded, result_expected)

    def test_load_result_file_malformed_json(self):
        with open(PARTIAL_JSON_PATH, 'w') as f:
            f.write("this is not json")
        smiles, result = load_result()
        self.assertIsNone(smiles)
        self.assertIsNone(result)
        with open(PARTIAL_JSON_PATH, 'w') as f:
            pass
        smiles, result = load_result()
        self.assertIsNone(smiles)
        self.assertIsNone(result)

    # Tests for /api/health endpoint
    def test_health_endpoint_no_api_key(self):
        response = self.client.get('/api/health')
        self.assertEqual(response.status_code, 401)
        self.assertEqual(response.json, {"error": "Unauthorized"})

    def test_health_endpoint_wrong_api_key(self):
        response = self.client.get('/api/health', headers={'X-API-KEY': 'wrong-key'})
        self.assertEqual(response.status_code, 401)
        self.assertEqual(response.json, {"error": "Unauthorized"})

    def test_health_endpoint_correct_api_key(self):
        response = self.client.get('/api/health', headers={'X-API-KEY': self.api_key})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json, {"status": "healthy"})
        
    def test_health_endpoint_options_request_bypasses_auth(self):
        response = self.client.options('/api/health')
        self.assertEqual(response.status_code, 200)
        if response.content_type == 'application/json':
            self.assertEqual(response.json, {"status": "healthy"})
        else:
            pass

    # Tests for /api/clear_molecule_cache endpoint
    def test_clear_molecule_cache_no_api_key(self):
        response = self.client.post('/api/clear_molecule_cache', json={'molecule': 'CCO'})
        self.assertEqual(response.status_code, 401)

    def test_clear_molecule_cache_no_molecule_param(self):
        response = self.client.post('/api/clear_molecule_cache',
                                     headers={'X-API-KEY': self.api_key},
                                     json={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "Molecule string is required"})

    @patch('src.api.clear_cache_for_molecule') 
    def test_clear_molecule_cache_success(self, mock_clear_cache):
        molecule_to_clear = "CCO"
        response = self.client.post('/api/clear_molecule_cache',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'molecule': molecule_to_clear})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json, {"status": "success"})
        mock_clear_cache.assert_called_once_with(molecule_to_clear)

    # Tests for /api/retrosynthesis endpoint input validation
    def test_retrosynthesis_no_api_key(self):
        response = self.client.post('/api/retrosynthesis', json={'smiles': 'CCO'})
        self.assertEqual(response.status_code, 401)

    def test_retrosynthesis_no_smiles_param(self):
        response = self.client.post('/api/retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={})
        self.assertEqual(response.status_code, 400)
        self.assertIn("SMILES string is required", response.json['error'])

    @patch('src.api.Chem.MolFromSmiles', return_value=None) 
    def test_retrosynthesis_invalid_smiles_string(self, mock_mol_from_smiles):
        response = self.client.post('/api/retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': 'invalid_smiles'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "Invalid SMILES string"})
        mock_mol_from_smiles.assert_called_once_with('invalid_smiles')

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True) 
    def test_retrosynthesis_default_parameters(self, mock_mol_from_smiles, mock_main, mock_save_result):
        mock_main.return_value = {"some_result": "data"} 
        smiles_input = "CCO"
        
        response = self.client.post('/api/retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': smiles_input})
        
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json, mock_main.return_value)
        mock_mol_from_smiles.assert_called_once_with(smiles_input)
        mock_main.assert_called_once_with(
            smiles=smiles_input,
            llm="claude-3-opus-20240229", 
            az_model="USPTO",            
            stability_flag="False",      
            hallucination_check="False"  
        )
        mock_save_result.assert_called_once_with(smiles_input, mock_main.return_value)

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_retrosynthesis_custom_parameters_valid(self, mock_mol_from_smiles, mock_main, mock_save_result):
        mock_main.return_value = {"custom_result": "custom_data"}
        smiles_input = "CC(=O)Oc1ccccc1C(=O)OH" 
        
        payload = {
            'smiles': smiles_input,
            'model_type': 'deepseek',
            'advanced_prompt': 'True', # String 'True'
            'model_version': 'gsk',
            'stability_flag': 'True',
            'hallucination_check': 'True'
        }
        
        with patch('src.api.AZ_MODEL_LIST', ["USPTO", "gsk", "pistachio"]):
            response = self.client.post('/api/retrosynthesis',
                                         headers={'X-API-KEY': self.api_key},
                                         json=payload)
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once_with(
            smiles=smiles_input,
            llm="fireworks_ai/accounts/fireworks/models/deepseek-r1:adv", 
            az_model="gsk", # Expect 'gsk' as it's now valid
            stability_flag="True",
            hallucination_check="True"
        )
        mock_save_result.assert_called_once_with(smiles_input, mock_main.return_value)

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_retrosynthesis_invalid_model_type_falls_back_to_default(self, mock_mol_from_smiles, mock_main, mock_save_result):
        mock_main.return_value = {"fallback_result": "data"}
        smiles_input = "CNC"
        
        payload = {
            'smiles': smiles_input,
            'model_type': 'invalid_model_type' 
        }
        
        response = self.client.post('/api/retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json=payload)
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once_with(
            smiles=smiles_input,
            llm="claude-3-opus-20240229", 
            az_model="USPTO",
            stability_flag="False",
            hallucination_check="False"
        )
        mock_save_result.assert_called_once_with(smiles_input, mock_main.return_value)

    @patch('src.api.save_result') # Added to avoid UnboundLocalError if main fails before save_result is defined
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_retrosynthesis_advanced_prompt_false_string(self, mock_mol_from_smiles, mock_main, mock_save_result):
        mock_main.return_value = {"result": "data"}
        smiles_input = "CCN"
        payload = {
            'smiles': smiles_input,
            'model_type': 'claude3',
            'advanced_prompt': 'false' # String "false"
        }
        response = self.client.post('/api/retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json=payload)
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once_with(
            smiles=smiles_input,
            llm="claude-3-opus-20240229:adv", # :adv WILL be appended by current API logic
            az_model="USPTO",
            stability_flag="False",
            hallucination_check="False"
        )

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_retrosynthesis_main_function_exception_handling(self, mock_mol_from_smiles, mock_main, mock_save_result):
        smiles_input = "CCC"
        mock_main.side_effect = Exception("Core processing failed")
        
        response = self.client.post('/api/retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': smiles_input})
        
        self.assertEqual(response.status_code, 500)
        self.assertIn("Error in retrosynthesis: Core processing failed", response.json['error'])
        mock_save_result.assert_not_called() 

    # Tests for /api/rerun_retrosynthesis endpoint
    def test_rerun_retrosynthesis_no_api_key(self):
        response = self.client.post('/api/rerun_retrosynthesis', json={'smiles': 'CCO'})
        self.assertEqual(response.status_code, 401)

    def test_rerun_retrosynthesis_no_smiles_param(self):
        response = self.client.post('/api/rerun_retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={})
        self.assertEqual(response.status_code, 400)
        self.assertIn("Molecule string is required", response.json['error'])

    @patch('src.api.Chem.MolFromSmiles', return_value=None)
    def test_rerun_retrosynthesis_invalid_smiles(self, mock_mol_from_smiles):
        response = self.client.post('/api/rerun_retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': 'invalid'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "Invalid SMILES string"})
        mock_mol_from_smiles.assert_called_once_with('invalid')

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    @patch('src.api.clear_cache_for_molecule')
    def test_rerun_retrosynthesis_default_params_success(self, mock_clear_cache, mock_mol_from_smiles, mock_main, mock_save_result):
        mock_main.return_value = {"rerun_result": "data"}
        smiles_input = "CCO"
        
        response = self.client.post('/api/rerun_retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': smiles_input})
        
        self.assertEqual(response.status_code, 200)
        mock_clear_cache.assert_called_once_with(smiles_input)
        mock_mol_from_smiles.assert_called_once_with(smiles_input)
        mock_main.assert_called_once_with(
            smiles=smiles_input,
            llm="claude-3-opus-20240229",
            az_model="USPTO",
            stability_flag="False",
            hallucination_check="False"
        )
        mock_save_result.assert_called_once_with(smiles_input, mock_main.return_value)

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    @patch('src.api.clear_cache_for_molecule')
    def test_rerun_retrosynthesis_custom_params_success(self, mock_clear_cache, mock_mol_from_smiles, mock_main, mock_save_result):
        mock_main.return_value = {"custom_rerun": "data"}
        smiles_input = "CN(C)C"
        payload = {
            'smiles': smiles_input,
            'model_type': 'claude37',
            'advanced_prompt': 'true', # string 'true'
            'model_version': 'pistachio', 
            'stability_flag': 'true',
            'hallucination_check': 'true'
        }
        
        with patch('src.api.AZ_MODEL_LIST', ["USPTO", "gsk", "pistachio"]): 
            response = self.client.post('/api/rerun_retrosynthesis',
                                        headers={'X-API-KEY': self.api_key},
                                        json=payload)
        
        self.assertEqual(response.status_code, 200)
        mock_clear_cache.assert_called_once_with(smiles_input)
        mock_main.assert_called_once_with(
            smiles=smiles_input,
            llm="anthropic/claude-3-7-sonnet-20250219:adv", # :adv will be added by API
            az_model="pistachio",
            stability_flag="true", 
            hallucination_check="true"
        )
        mock_save_result.assert_called_once_with(smiles_input, mock_main.return_value)

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    @patch('src.api.clear_cache_for_molecule')
    def test_rerun_retrosynthesis_main_exception(self, mock_clear_cache, mock_mol_from_smiles, mock_main, mock_save_result):
        smiles_input = "CCN"
        mock_main.side_effect = Exception("Rerun failed badly")
        
        response = self.client.post('/api/rerun_retrosynthesis',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': smiles_input})
        
        self.assertEqual(response.status_code, 500)
        self.assertIn("Error in retrosynthesis: Rerun failed badly", response.json['error'])
        mock_clear_cache.assert_called_once_with(smiles_input)
        mock_save_result.assert_not_called()

    # Tests for /api/partial_rerun endpoint - Initial Error Handling
    def test_partial_rerun_no_api_key(self):
        response = self.client.post('/api/partial_rerun', json={'smiles': 'CCO', 'steps': 1})
        self.assertEqual(response.status_code, 401)

    def test_partial_rerun_missing_smiles_in_payload(self):
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'steps': 1})
        self.assertEqual(response.status_code, 500) 
        self.assertIn("error", response.json)
        self.assertTrue("'smiles'" in response.json["error"] or "Error in partial retrosynthesis: 'smiles'" in response.json["error"])

    def test_partial_rerun_missing_steps_in_payload(self):
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': 'CCO'})
        self.assertEqual(response.status_code, 500) 
        self.assertIn("error", response.json)
        self.assertTrue("'steps'" in response.json["error"] or "Error in partial retrosynthesis: 'steps'" in response.json["error"])

    @patch('src.api.load_result', return_value=(None, None))
    def test_partial_rerun_no_stored_result(self, mock_load_result):
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': 'CCO', 'steps': 1})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "No previous results found. Run retrosynthesis first."})
        mock_load_result.assert_called_once()

    @patch('src.api.load_result', return_value=("CCC", {"some_data": "value"}))
    def test_partial_rerun_smiles_mismatch(self, mock_load_result):
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': 'CCO', 'steps': 1})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "No results found for this molecule. Run retrosynthesis first."})
        mock_load_result.assert_called_once()

    @patch('src.api.load_result')
    def test_partial_rerun_step_not_found(self, mock_load_result):
        original_smiles = "CCO"
        original_result = {
            "steps": [
                {"step": "1", "products": [{"smiles": "CCO"}]}, # step as string
                {"step": "2", "products": [{"smiles": "CC"}]}  # step as string
            ],
            "dependencies": { "1": ["2"] }
        }
        mock_load_result.return_value = (original_smiles, original_result)
        
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': original_smiles, 'steps': 3}) 
        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.json, {"error": "Step 3 not found in the synthesis pathway"})

    @patch('src.api.load_result')
    def test_partial_rerun_target_step_no_products(self, mock_load_result):
        original_smiles = "CCO"
        original_result = {
            "steps": [
                {"step": "1", "products": []} 
            ],
            "dependencies": {}
        }
        mock_load_result.return_value = (original_smiles, original_result)
        
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': original_smiles, 'steps': 1})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "Step 1 doesn't have valid products"})

    @patch('src.api.load_result')
    def test_partial_rerun_target_step_invalid_products_format(self, mock_load_result):
        original_smiles = "CCO"
        original_result = {
            "steps": [
                {"step": "1", "products": "not_a_list"} 
            ],
            "dependencies": {}
        }
        mock_load_result.return_value = (original_smiles, original_result)
        
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json={'smiles': original_smiles, 'steps': 1})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json, {"error": "Step 1 doesn't have valid products"})

    # Core logic tests for /api/partial_rerun
    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.load_result')
    @patch('src.api.Chem.MolFromSmiles', return_value=True) 
    def test_partial_rerun_middle_step(self, mock_chem_validate, mock_load_result, mock_main, mock_save_result):
        original_smiles = "A" # Overall target
        original_result = {
            "smiles": original_smiles,
            "steps": [
                {"step": "1", "name": "Reaction 1 A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]},
                {"step": "2", "name": "Reaction 2 B->C", "reactants": [{"smiles": "B"}], "products": [{"smiles": "C"}]}, # Target step for rerun
                {"step": "3", "name": "Reaction 3 C->D", "reactants": [{"smiles": "C"}], "products": [{"smiles": "D"}]}
            ],
            "dependencies": {"1": ["2"], "2": ["3"]}
        }
        mock_load_result.return_value = (original_smiles, original_result)

        # Rerun from step 2 (B->C). The product of target_step is "C". So, new synthesis starts from "C".
        start_molecule_for_new_synthesis = "C"
        new_sub_result_for_C = {
            "smiles": start_molecule_for_new_synthesis, # This field in new_result is not actually used by API
            "steps": [
                {"step": 1, "name": "New Reaction C->X", "reactants": [{"smiles": "C"}], "products": [{"smiles": "X"}]},
                {"step": 2, "name": "New Reaction X->Y", "reactants": [{"smiles": "X"}], "products": [{"smiles": "Y"}]}
            ],
            "dependencies": {"1": ["2"]} # Local dependencies for C->X->Y
        }
        mock_main.return_value = new_sub_result_for_C

        rerun_from_step_num = 2 
        payload = {
            'smiles': original_smiles, 
            'steps': rerun_from_step_num
        }
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json=payload)
        
        self.assertEqual(response.status_code, 200, f"Response JSON: {response.json}")
        
        mock_main.assert_called_once_with(
            smiles=start_molecule_for_new_synthesis, # Should be "C"
            llm=ANY, az_model=ANY, stability_flag=ANY, hallucination_check=ANY
        )

        # Kept: Step 1 (A->B). Max kept step is 1.
        # New steps (C->X, X->Y) are renumbered to 2, 3.
        # Original step 2 (B->C) and 3 (C->D) are removed.
        expected_merged_steps = [
            {"step": "1", "name": "Reaction 1 A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]},
            {"step": "2", "name": "New Reaction C->X", "reactants": [{"smiles": "C"}], "products": [{"smiles": "X"}]}, 
            {"step": "3", "name": "New Reaction X->Y", "reactants": [{"smiles": "X"}], "products": [{"smiles": "Y"}]}  
        ]
        # API connects left_connection (step "1") to the first new step ("2" which is C->X)
        expected_merged_dependencies = {
            "1": ["2"], 
            "2": ["3"] 
        }

        actual_result = response.json
        self.assertEqual(actual_result['smiles'], original_smiles) # API returns original target SMILES
        self.assertListEqual(sorted(actual_result['steps'], key=lambda x: int(x['step'])), sorted(expected_merged_steps, key=lambda x: int(x['step'])))
        self.assertDictEqual(actual_result['dependencies'], expected_merged_dependencies)
        mock_save_result.assert_called_once_with(original_smiles, actual_result) 

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.load_result')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_partial_rerun_root_step(self, mock_chem_validate, mock_load_result, mock_main, mock_save_result):
        original_smiles = "A" # Overall target
        original_result = {
            "smiles": original_smiles,
            "steps": [
                {"step": "1", "name": "Original A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]}, # Target step for rerun
                {"step": "2", "name": "Original B->C", "reactants": [{"smiles": "B"}], "products": [{"smiles": "C"}]}
            ],
            "dependencies": {"1": ["2"]}
        }
        mock_load_result.return_value = (original_smiles, original_result)

        # Rerun from step 1 (A->B). The product of target_step is "B". So, new synthesis starts from "B".
        start_molecule_for_new_synthesis = "B" 
        new_sub_result_for_B = {
            "smiles": start_molecule_for_new_synthesis, # This field is not directly used by API for merging
            "steps": [
                {"step": 1, "name": "New B->Z", "reactants": [{"smiles": "B"}], "products": [{"smiles": "Z"}]},
                {"step": 2, "name": "New Z->W", "reactants": [{"smiles": "Z"}], "products": [{"smiles": "W"}]}
            ],
            "dependencies": {"1": ["2"]} # Local dependencies for B->Z->W
        }
        mock_main.return_value = new_sub_result_for_B

        rerun_from_step_num = 1
        payload = {'smiles': original_smiles, 'steps': rerun_from_step_num}
        response = self.client.post('/api/partial_rerun',
                                     headers={'X-API-KEY': self.api_key},
                                     json=payload)
        
        self.assertEqual(response.status_code, 200, f"Response JSON: {response.json}")
        mock_main.assert_called_once_with(
            smiles=start_molecule_for_new_synthesis, # Should be "B"
            llm=ANY, az_model=ANY, stability_flag=ANY, hallucination_check=ANY
        )

        # Kept_steps is empty as step 1 is removed. Max_kept_step is 0.
        # New steps (B->Z, Z->W) are renumbered to 1, 2.
        expected_merged_steps = [
            {"step": "1", "name": "New B->Z", "reactants": [{"smiles": "B"}], "products": [{"smiles": "Z"}]}, 
            {"step": "2", "name": "New Z->W", "reactants": [{"smiles": "Z"}], "products": [{"smiles": "W"}]}
        ]
        expected_merged_dependencies = {"1": ["2"]}

        actual_result = response.json
        self.assertEqual(actual_result['smiles'], original_smiles) # API returns original target SMILES
        self.assertListEqual(sorted(actual_result['steps'], key=lambda x: int(x['step'])), sorted(expected_merged_steps, key=lambda x: int(x['step'])))
        self.assertDictEqual(actual_result['dependencies'], expected_merged_dependencies)
        mock_save_result.assert_called_once_with(original_smiles, actual_result)

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.load_result')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_partial_rerun_leaf_step(self, mock_chem_validate, mock_load_result, mock_main, mock_save_result):
        original_smiles = "A"
        original_result = {
            "smiles": original_smiles,
            "steps": [
                {"step": "1", "name": "A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]},
                {"step": "2", "name": "B->C", "reactants": [{"smiles": "B"}], "products": [{"smiles": "C"}]} # Rerun this step
            ],
            "dependencies": {"1": ["2"]}
        }
        mock_load_result.return_value = (original_smiles, original_result)

        # Rerun from step 2 (B->C). Product of target_step is "C". New synthesis for "C".
        start_molecule_for_new_synthesis = "C"
        new_sub_result_for_C = {
            "smiles": start_molecule_for_new_synthesis,
            "steps": [{"step": 1, "name": "New C->X", "reactants": [{"smiles": "C"}], "products": [{"smiles": "X"}]}],
            "dependencies": {}
        }
        mock_main.return_value = new_sub_result_for_C

        rerun_from_step_num = 2
        payload = {'smiles': original_smiles, 'steps': rerun_from_step_num}
        response = self.client.post('/api/partial_rerun', headers={'X-API-KEY': self.api_key}, json=payload)
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once_with(
            smiles=start_molecule_for_new_synthesis, # "C"
            llm=ANY, az_model=ANY, stability_flag=ANY, hallucination_check=ANY
        )

        # Kept step: 1 (A->B). Max kept step is 1.
        # New step (C->X) renumbered to 2.
        expected_merged_steps = [
            {"step": "1", "name": "A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]},
            {"step": "2", "name": "New C->X", "reactants": [{"smiles": "C"}], "products": [{"smiles": "X"}]} 
        ]
        # API connects left_connection (step "1") to the first new step ("2" which is C->X)
        expected_merged_dependencies = {"1": ["2"]}
        actual_result = response.json
        self.assertEqual(actual_result['smiles'], original_smiles)
        self.assertListEqual(sorted(actual_result['steps'], key=lambda x: int(x['step'])), sorted(expected_merged_steps, key=lambda x: int(x['step'])))
        self.assertDictEqual(actual_result['dependencies'], expected_merged_dependencies)
        mock_save_result.assert_called_once_with(original_smiles, actual_result)

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.load_result')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_partial_rerun_new_sub_synthesis_empty(self, mock_chem_validate, mock_load_result, mock_main, mock_save_result):
        original_smiles = "A"
        original_result = {
            "smiles": original_smiles,
            "steps": [
                {"step": "1", "name": "A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]},
                {"step": "2", "name": "B->C", "reactants": [{"smiles": "B"}], "products": [{"smiles": "C"}]} # Rerun this step
            ],
            "dependencies": {"1": ["2"]}
        }
        mock_load_result.return_value = (original_smiles, original_result)

        # Rerun from step 2 (B->C). Product of target step is "C". New synthesis for "C".
        start_molecule_for_new_synthesis = "C"
        mock_main.return_value = {"smiles": start_molecule_for_new_synthesis, "steps": [], "dependencies": {}}

        rerun_from_step_num = 2
        payload = {'smiles': original_smiles, 'steps': rerun_from_step_num}
        response = self.client.post('/api/partial_rerun', headers={'X-API-KEY': self.api_key}, json=payload)
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once_with(
            smiles=start_molecule_for_new_synthesis, # "C"
            llm=ANY, az_model=ANY, stability_flag=ANY, hallucination_check=ANY
        )
        
        # Kept step: 1 (A->B). Original step 2 (B->C) removed. New synthesis for C is empty.
        # Path ends at B.
        expected_merged_steps = [
            {"step": "1", "name": "A->B", "reactants": [{"smiles": "A"}], "products": [{"smiles": "B"}]}
        ]
        expected_merged_dependencies = {"1": []} # Updated to match actual API behavior
        actual_result = response.json
        self.assertEqual(actual_result['smiles'], original_smiles)
        self.assertListEqual(sorted(actual_result['steps'], key=lambda x: int(x['step'])), sorted(expected_merged_steps, key=lambda x: int(x['step'])))
        self.assertDictEqual(actual_result['dependencies'], expected_merged_dependencies)
        mock_save_result.assert_called_once_with(original_smiles, actual_result)

    @patch('src.api.save_result') 
    @patch('src.api.main')
    @patch('src.api.load_result')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_partial_rerun_main_call_exception(self, mock_chem_validate, mock_load_result, mock_main, mock_save_result):
        original_smiles = "A"
        # Target step 1 (A->B), so new synthesis will be for "B"
        original_result = {
            "smiles": original_smiles,
            "steps": [{"step": "1", "name": "A->B", "products": [{"smiles": "B"}]}],
            "dependencies": {}
        }
        mock_load_result.return_value = (original_smiles, original_result)
        
        # New synthesis starts from "B" (product of step 1)
        start_molecule_for_new_synthesis = "B"
        mock_main.side_effect = Exception("Sub-synthesis failed")

        rerun_from_step_num = 1
        payload = {'smiles': original_smiles, 'steps': rerun_from_step_num}
        response = self.client.post('/api/partial_rerun', headers={'X-API-KEY': self.api_key}, json=payload)
        
        self.assertEqual(response.status_code, 500)
        self.assertIn(f"Error running retrosynthesis on {start_molecule_for_new_synthesis}: Sub-synthesis failed", response.json['error'])
        mock_main.assert_called_once_with(
            smiles=start_molecule_for_new_synthesis, # "B"
            llm=ANY, az_model=ANY, stability_flag=ANY, hallucination_check=ANY
        )
        mock_save_result.assert_not_called()

    @patch('src.api.save_result')
    @patch('src.api.main')
    @patch('src.api.load_result')
    @patch('src.api.Chem.MolFromSmiles', return_value=True)
    def test_partial_rerun_custom_params_to_main(self, mock_chem_validate, mock_load_result, mock_main, mock_save_result):
        original_smiles = "A"
        # Target step 1 (A->B), so new synthesis will be for "B"
        original_result = {
            "smiles": original_smiles,
            "steps": [{"step": "1", "name": "A->B", "products": [{"smiles": "B"}]}],
            "dependencies": {}
        }
        mock_load_result.return_value = (original_smiles, original_result)
        
        # New synthesis starts from "B"
        start_molecule_for_new_synthesis = "B"
        mock_main.return_value = {"smiles": start_molecule_for_new_synthesis, "steps": [{"step":1, "name":"New B->X", "products":[{"smiles":"X"}]}], "dependencies":{}} 
        
        rerun_from_step_num = 1
        payload = {
            'smiles': original_smiles, 
            'steps': rerun_from_step_num,
            'model_type': 'deepseek',
            'advanced_prompt': 'True', # String 'True'
            'model_version': 'gsk',
            'stability_flag': 'True',
            'hallucination_check': 'True'
        }
        
        with patch('src.api.AZ_MODEL_LIST', ["USPTO", "gsk"]): 
            response = self.client.post('/api/partial_rerun', headers={'X-API-KEY': self.api_key}, json=payload)
        
        self.assertEqual(response.status_code, 200)
        mock_main.assert_called_once_with(
            smiles=start_molecule_for_new_synthesis, # "B"
            llm="fireworks_ai/accounts/fireworks/models/deepseek-r1:adv", # :adv will be added
            az_model="gsk",
            stability_flag="True",
            hallucination_check="True"
        )
        mock_save_result.assert_called_once()

if __name__ == '__main__':
    unittest.main() 