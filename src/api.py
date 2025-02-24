from flask import Flask, request, jsonify
from flask_cors import CORS
from functools import wraps
from dotenv import load_dotenv
from rdkit import Chem

import rootutils

import json                 # For handling JSON data
import traceback           # For detailed error tracking
root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

load_dotenv()
from src.main import main
from src.cache import clear_cache_for_molecule
from src.variables import AZ_MODEL_LIST

app = Flask(__name__)
CORS(app)

# Predefined API key for authentication
API_KEY = "your-secure-api-key"

# Global storage for the latest retrosynthesis results
latest_results = {}

# Authentication decorator
def require_api_key(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        # Allow OPTIONS requests to bypass API key check
        if request.method == 'OPTIONS':
            return f(*args, **kwargs)
            
        api_key = request.headers.get('X-API-KEY')
        if api_key and api_key == API_KEY:
            return f(*args, **kwargs)
        else:
            return jsonify({"error": "Unauthorized"}), 401
    return decorated_function

@app.route('/api/retrosynthesis', methods=['POST'])
@require_api_key
def retrosynthesis_api():
    """
    Endpoint to perform retrosynthesis on a SMILES string.
    """
    global latest_results
    
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({
            "error":
            "SMILES string is required. Please include a 'smiles' field"
        }), 400

    smiles = data['smiles']

    # Check if the SMILES string is valid
    if not Chem.MolFromSmiles(smiles):
        return jsonify({"error": "Invalid SMILES string"}), 400

    # -----------------
    # Advanced model - DeepSeek-R1
    deepseek_r1 = False
    try:
        advanced_model: str = data['advanced_model']
        if advanced_model.lower() == "true":
            deepseek_r1 = True
    except Exception as e:
        print(e)
        advanced_model = False

    if deepseek_r1:
        llm = "deepinfra/deepseek-ai/DeepSeek-R1"
    else:
        llm = "claude-3-opus-20240229"

    # -----------------
    # Advanced Prompt - To use the more guardrails prompt
    advanced_prompt = False
    try:
        advanced_prompt: str = data['advanced_prompt']
        if advanced_prompt.lower() == "true":
            advanced_prompt = True
    except Exception as e:
        print(e)
        advanced_prompt = False

    if advanced_prompt:
        llm = llm + ":adv"

    # -----------------
    # Choose AiZynthFinder model
    az_model = "USPTO"
    try:
        az_model: str = data['model_version']
        assert az_model in AZ_MODEL_LIST
    except Exception as e:
        print(e)
        az_model = "USPTO"

    # -----------------
    # Run retrosynthesis
    try:
        result = main(smiles=smiles, llm=llm, az_model=az_model)
        # Store the result for potential partial reruns
        latest_results[smiles] = result
    except Exception as e:
        print(e)
        return jsonify({"error": "Error in retrosynthesis, Please rerun"}), 500
    return jsonify(result), 200


@app.route('/api/health', methods=['GET'])
@require_api_key
def health():
    """
    Endpoint to check the health of the API.
    """
    return jsonify({"status": "healthy"}), 200


@app.route('/api/clear_molecule_cache', methods=['POST'])
@require_api_key
def clear_molecule_cache():
    """
    Endpoint to clear the cache for a specific molecule.
    """
    data = request.get_json()
    if not data or 'molecule' not in data:
        return jsonify({"error": "Molecule string is required"}), 400

    molecule = data['molecule']
    clear_cache_for_molecule(molecule)
    return jsonify({"status": "success"}), 200


@app.route('/api/rerun_retrosynthesis', methods=['POST'])
@require_api_key
def rerun_retrosynthesis():
    """
    Endpoint to rerun retrosynthesis for a specific molecule.
    """
    global latest_results
    
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({
            "error":
            "Molecule string is required, Please include a 'smiles' field"
        }), 400

    molecule = data['smiles']

    # Clear the cache for the molecule
    clear_cache_for_molecule(molecule)
    deepseek_r1 = False
    try:
        advanced_model: str = data['advanced_model']
        if advanced_model.lower() == "true":
            deepseek_r1 = True
    except Exception as e:
        print(e)
        advanced_model = False

    if not Chem.MolFromSmiles(molecule):
        return jsonify({"error": "Invalid SMILES string"}), 400

    if deepseek_r1:
        llm = "deepinfra/deepseek-ai/DeepSeek-R1"
    else:
        llm = "claude-3-opus-20240229"

    # Advanced prompt handling
    advanced_prompt = False
    try:
        advanced_prompt: str = data['advanced_prompt']
        if advanced_prompt.lower() == "true":
            advanced_prompt = True
    except Exception as e:
        print(e)
        advanced_prompt = False

    if advanced_prompt:
        llm = llm + ":adv"

    # Choose AiZynthFinder model
    az_model = "USPTO"
    try:
        az_model: str = data['model_version']
        assert az_model in AZ_MODEL_LIST
    except Exception as e:
        print(e)
        az_model = "USPTO"

    # Rerun retrosynthesis
    try:
        result = main(smiles=molecule, llm=llm, az_model=az_model)
        # Store the result for potential partial reruns
        latest_results[molecule] = result
    except Exception as e:
        print(e)
        return jsonify({"error": "Error in retrosynthesis, Please rerun"}), 500
    return jsonify(result), 200

@app.route('/api/partial_rerun', methods=['POST'])
@require_api_key
def partial_rerun():
    """
    Endpoint to partially rerun retrosynthesis from a specific step.
    Uses the stored results from the most recent retrosynthesis run.
    """
    global latest_results
    
    print("\n=== Starting Partial Rerun Process ===")
    
    data = request.get_json()
    print(f"Received request data: {json.dumps(data, indent=2)}")
    
    try:
        smiles = data['smiles']
        from_step = int(data['steps'])
        
        # Check if we have stored results for this molecule
        if smiles not in latest_results:
            return jsonify({"error": "No previous results found for this molecule. Run retrosynthesis first."}), 400
            
        # Get the original result from our stored data
        original_result = latest_results[smiles]
        print(f"Found stored result for SMILES: {smiles}")
        
        # Get the starting molecule from the specified step
        target_step = next(
            (step for step in original_result['steps'] 
             if int(step['step']) == from_step),
            None
        )
        
        if not target_step:
            return jsonify({"error": f"Step {from_step} not found in the synthesis pathway"}), 404
            
        start_molecule = target_step['reactants'][0]['smiles']
        print(f"\nStarting new synthesis from molecule: {start_molecule}")
        
        # Run new synthesis on the starting molecule by calling the retrosynthesis endpoint
        retro_request = {
            'smiles': start_molecule,
            'advanced_model': data.get('advanced_model', "false"),
            'advanced_prompt': data.get('advanced_prompt', "false"),
            'model_version': data.get('model_version', "USPTO")
        }
        
        with app.test_client() as client:
            new_result_response = client.post(
                '/api/retrosynthesis',
                json=retro_request,
                headers={'X-API-KEY': API_KEY}
            )
            if new_result_response.status_code != 200:
                return jsonify({"error": f"Error running retrosynthesis on molecule {start_molecule}"}), 500
            new_result = new_result_response.get_json()
        
        # Analyze the original structure to determine which steps to keep
        # First identify all steps that are downstream of the target step
        downstream_steps = set()
        
        def find_downstream_steps(step_num, deps):
            """Recursively identify all steps that depend on the given step"""
            for s, dep_list in deps.items():
                if str(step_num) in dep_list:
                    downstream_steps.add(s)
                    find_downstream_steps(s, deps)
        
        find_downstream_steps(from_step, original_result['dependencies'])
        
        # Identify all steps in the branch containing the target step
        target_branch = {str(from_step)}
        
        def find_branch_steps(step_num, deps):
            """Recursively identify all steps in the same branch"""
            if step_num in deps:
                for dep in deps[step_num]:
                    target_branch.add(dep)
                    find_branch_steps(dep, deps)
        
        find_branch_steps(str(from_step), original_result['dependencies'])
        
        # Steps to remove = steps in target branch + downstream steps
        steps_to_remove = target_branch.union(downstream_steps)
        
        # Keep steps that are not in the steps_to_remove set
        kept_steps = [
            step.copy() 
            for step in original_result['steps'] 
            if str(step['step']) not in steps_to_remove
        ]
        
        # Keep dependencies for steps we're keeping, removing any references to removed steps
        kept_deps = {}
        for step_num, deps in original_result['dependencies'].items():
            if step_num not in steps_to_remove:
                # Filter out dependencies that are in steps_to_remove
                kept_deps[step_num] = [d for d in deps if d not in steps_to_remove]
        
        # Find max step number from kept steps
        max_step = 0
        if kept_steps:
            max_step = max(int(step['step']) for step in kept_steps)
        
        # Find the parent step to which we need to attach the new branch
        parent_step = None
        for step_num, deps in original_result['dependencies'].items():
            if str(from_step) in deps and step_num not in steps_to_remove:
                parent_step = step_num
                break
        
        # Adjust new steps (renumber them starting from max_step + 1)
        new_steps = []
        step_mapping = {}
        
        for idx, step in enumerate(new_result['steps']):
            new_step_num = max_step + 1 + idx
            step_mapping[step['step']] = str(new_step_num)
            
            adjusted_step = step.copy()
            adjusted_step['step'] = str(new_step_num)
            new_steps.append(adjusted_step)
        
        # Adjust dependencies for new steps
        new_deps = {}
        for old_num, deps in new_result['dependencies'].items():
            new_num = step_mapping[old_num]
            new_deps[new_num] = [step_mapping[d] for d in deps]
        
        # Link first new step to parent step if one exists
        if parent_step:
            if parent_step in kept_deps:
                kept_deps[parent_step].append(new_steps[0]['step'])
            else:
                kept_deps[parent_step] = [new_steps[0]['step']]
        
        # Merge results
        merged_result = {
            'steps': kept_steps + new_steps,
            'dependencies': {**kept_deps, **new_deps}
        }
        
        # Store the merged result as the latest result for this molecule
        latest_results[smiles] = merged_result
        
        print("\n=== Partial Rerun Complete ===")
        return jsonify(merged_result), 200
        
    except Exception as e:
        print(f"\nERROR in partial rerun:")
        print(f"Exception: {str(e)}")
        print(f"Traceback: {traceback.format_exc()}")
        return jsonify({
            "error": f"Error in partial retrosynthesis: {str(e)}"
        }), 500

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)