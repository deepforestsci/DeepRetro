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
    except Exception as e:
        print(e)
        return jsonify({"error": "Error in retrosynthesis, Please rerun"}), 500
    return jsonify(result), 200

@app.route('/api/partial_rerun', methods=['POST'])
@require_api_key
def partial_rerun():
   print("\n=== Starting Partial Rerun Process ===")
   
   data = request.get_json()
   print(f"Received request data: {json.dumps(data, indent=2)}")
   
   try:
       smiles = data['smiles']
       from_step = int(data['steps'])
       pathway_id = data.get('pathway_id', 1)
       
       # Get original result via retrosynthesis endpoint
       retro_request = {
           'smiles': smiles,
           'advanced_model': data.get('advanced_model', "false"),
           'advanced_prompt': data.get('advanced_prompt', "false"),
           'model_version': data.get('model_version', "USPTO")
       }
       
       with app.test_client() as client:
           original_response = client.post(
               '/api/retrosynthesis',
               json=retro_request,
               headers={'X-API-KEY': API_KEY}
           )
           if original_response.status_code != 200:
               return jsonify({"error": "Could not get original synthesis"}), 404
           original_result = original_response.get_json()
       
       # Get starting molecule
       target_step = next(
           (step for step in original_result['steps'] 
            if int(step['step']) == from_step),
           None
       )
       if not target_step:
           return jsonify({"error": f"Step {from_step} not found"}), 404
           
       start_molecule = target_step['reactants'][0]['smiles']
       print(f"\nStarting new synthesis from molecule: {start_molecule}")
       
       # Run new synthesis
       llm = ("deepinfra/deepseek-ai/DeepSeek-R1" 
              if data.get('advanced_model', "false").lower() == "true" 
              else "claude-3-opus-20240229")
       if data.get('advanced_prompt', "false").lower() == "true":
           llm += ":adv"

       new_result = main(
           smiles=start_molecule,
           llm=llm,
           az_model=data.get('model_version', "USPTO")
       )
       
       if isinstance(new_result, tuple):
           new_result = {'steps': new_result[0], 'dependencies': new_result[1]}
       
       # Keep steps before from_step
       kept_steps = [
           step.copy() 
           for step in original_result['steps'] 
           if int(step['step']) < from_step
       ]
       
       # Keep only dependencies for steps we're keeping
       kept_deps = {
           step_num: []  # Initialize with empty list - we'll set correct dependencies later
           for step_num in [str(i) for i in range(1, from_step)]
       }
       
       # Find max step number from kept steps
       max_step = max(int(step['step']) for step in kept_steps)
       
       # Adjust new steps
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
       
       # Link first new step to parent
       parent_step = str(from_step - 1)
       if parent_step in kept_deps:
           kept_deps[parent_step] = [new_steps[0]['step']]  # Set new step as only dependency
           
       # For all other kept steps, preserve original dependencies if they're within kept steps
       for step_num in kept_deps:
           if step_num != parent_step:  # Skip the parent step as we've already set its deps
               original_deps = original_result['dependencies'].get(step_num, [])
               kept_deps[step_num] = [d for d in original_deps if int(d) < from_step]
       
       # Merge results
       merged_result = {
           'steps': kept_steps + new_steps,
           'dependencies': {**kept_deps, **new_deps}
       }
       
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