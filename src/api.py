from flask import Flask, request, jsonify
from flask_cors import CORS
from functools import wraps
from dotenv import load_dotenv
from rdkit import Chem

import rootutils
import os
import hashlib
import json
import traceback

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


# Functions for JSON file storage
def save_result(smiles, result):
    """Save retrosynthesis result to a JSON file."""
    # Create results directory if it doesn't exist
    os.makedirs("results", exist_ok=True)

    # Create a safe filename using a hash of the SMILES
    filename = f"results/{hashlib.md5(smiles.encode()).hexdigest()}.json"

    # Save the result to the file
    with open(filename, 'w') as f:
        json.dump(result, f)

    print(f"Result saved to {filename}")
    return filename


def load_result(smiles):
    """Load retrosynthesis result from a JSON file."""
    filename = f"results/{hashlib.md5(smiles.encode()).hexdigest()}.json"

    # Check if the file exists
    if not os.path.exists(filename):
        print(f"No result file found for {smiles}")
        return None

    # Load the result from the file
    try:
        with open(filename, 'r') as f:
            result = json.load(f)
        print(f"Result loaded from {filename}")
        return result
    except Exception as e:
        print(f"Error loading result file: {e}")
        return None


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
        llm = "fireworks_ai/accounts/fireworks/models/deepseek-r1"
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
    # Stability check flag
    stability_flag = "False"
    try:
        stability_flag: str = data['stability_flag']
        assert stability_flag.lower() in ["false", "true"]
    except Exception as e:
        print(e)
        stability_flag = "False"

    # -----------------
    # Hallucination check flag
    hallucination_check = "False"
    try:
        hallucination_check: str = data['hallucination_check']
        assert hallucination_check.lower() in ["false", "true"]
    except Exception as e:
        print(e)
        hallucination_check = "False"

    # -----------------
    # Run retrosynthesis
    try:
        result = main(smiles=smiles,
                      llm=llm,
                      az_model=az_model,
                      stability_flag=stability_flag,
                      hallucination_check=hallucination_check)

        # Store the result in a JSON file
        save_result(smiles, result)

    except Exception as e:
        print(e)
        return jsonify(
            {"error":
             f"Error in retrosynthesis: {str(e)}. Please rerun."}), 500

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

    if not Chem.MolFromSmiles(molecule):
        return jsonify({"error": "Invalid SMILES string"}), 400

    # Clear the cache for the molecule
    clear_cache_for_molecule(molecule)

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
        llm = "fireworks_ai/accounts/fireworks/models/deepseek-r1"
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
    # Stability check flag
    stability_flag = "False"
    try:
        stability_flag: str = data['stability_flag']
        assert stability_flag.lower() in ["false", "true"]
    except Exception as e:
        print(e)
        stability_flag = "False"

    # -----------------
    # Hallucination check flag
    hallucination_check = "False"
    try:
        hallucination_check: str = data['hallucination_check']
        assert hallucination_check.lower() in ["false", "true"]
    except Exception as e:
        print(e)
        hallucination_check = "False"

    # Rerun retrosynthesis
    try:
        result = main(smiles=molecule,
                      llm=llm,
                      az_model=az_model,
                      stability_flag=stability_flag,
                      hallucination_check=hallucination_check)
        save_result(molecule, result)

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
    
    When rerunning a step, we remove that step and everything to its right in the synthesis pathway.
    """

    print("\n=== Starting Partial Rerun Process ===")

    data = request.get_json()
    print(f"Received request data: {json.dumps(data, indent=2)}")

    try:
        smiles = data['smiles']
        from_step = int(data['steps'])

        # Load previous results from JSON file
        original_result = load_result(smiles)
        if not original_result:
            return jsonify({
                "error":
                "No previous results found for this molecule. Run retrosynthesis first."
            }), 400

        print(f"Found stored result for SMILES: {smiles}")

        # Print the original dependency structure for debugging
        print(
            f"Original dependencies structure: {json.dumps(original_result.get('dependencies', {}), indent=2)}"
        )
        print(
            f"Original steps: {json.dumps([s['step'] for s in original_result.get('steps', [])])}"
        )

        # Get the starting molecule from the specified step
        target_step = next((step for step in original_result['steps']
                            if int(step['step']) == from_step), None)
        print(
            f"Original steps: {json.dumps([s['step'] for s in original_result['steps']], indent=2)}"
        )

        # Get the starting molecule from the specified step
        target_step = next((step for step in original_result['steps']
                            if int(step['step']) == from_step), None)

        if not target_step:
            return jsonify({
                "error":
                f"Step {from_step} not found in the synthesis pathway"
            }), 404

        start_molecule = target_step['reactants'][0]['smiles']
        print(f"\nStarting new synthesis from molecule: {start_molecule}")

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
            llm = "fireworks_ai/accounts/fireworks/models/deepseek-r1"
        else:
            llm = "claude-3-opus-20240229"

        # -----------------
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
        # Stability check flag
        stability_flag = "False"
        try:
            stability_flag: str = data['stability_flag']
            assert stability_flag.lower() in ["false", "true"]
        except Exception as e:
            print(e)
            stability_flag = "False"

        # -----------------
        # Hallucination check flag
        hallucination_check = "False"
        try:
            hallucination_check: str = data['hallucination_check']
            assert hallucination_check.lower() in ["false", "true"]
        except Exception as e:
            print(e)
            hallucination_check = "False"

        # Run new synthesis on the starting molecule
        try:
            new_result = main(smiles=start_molecule,
                              llm=llm,
                              az_model=az_model,
                              stability_flag=stability_flag,
                              hallucination_check=hallucination_check)
            print(
                f"New retrosynthesis result: {json.dumps(new_result, indent=2)}"
            )
        except Exception as e:
            print(
                f"Error running retrosynthesis on molecule {start_molecule}: {str(e)}"
            )
            return jsonify({
                "error":
                f"Error running retrosynthesis on {start_molecule}: {str(e)}"
            }), 500

        # The steps to remove are the target step and everything it depends on
        # In other words, the target step and everything to its right in the synthesis pathway
        steps_to_remove = {str(from_step)}
        steps_to_check = [str(from_step)]

        # Identify all steps that the target step depends on (to the right in the pathway)
        # This is the set of steps that will be replaced by the new synthesis
        while steps_to_check:
            current_step = steps_to_check.pop(0)
            if current_step in original_result.get('dependencies', {}):
                for dep_step in original_result['dependencies'][current_step]:
                    if dep_step not in steps_to_remove:
                        steps_to_remove.add(dep_step)
                        steps_to_check.append(dep_step)

        print(
            f"Steps to remove (target and everything to its right): {steps_to_remove}"
        )

        # We need to identify what step the target step is connected to on its left
        # This is where we'll connect the new synthesis
        left_connection = None
        for step, deps in original_result.get('dependencies', {}).items():
            if str(from_step) in deps and step not in steps_to_remove:
                left_connection = step
                break

        print(
            f"Left connection (step that the target connects to on the left): {left_connection}"
        )

        # Keep steps that are not in the steps_to_remove set
        kept_steps = [
            step.copy() for step in original_result['steps']
            if str(step['step']) not in steps_to_remove
        ]

        # Keep dependencies for steps we're keeping, removing any references to removed steps
        kept_deps = {}
        for step_num, deps in original_result.get('dependencies', {}).items():
            if step_num not in steps_to_remove:
                # Filter out dependencies that are in steps_to_remove
                kept_deps[step_num] = [
                    d for d in deps if d not in steps_to_remove
                ]

        print(f"Kept steps: {[s['step'] for s in kept_steps]}")
        print(f"Kept dependencies: {kept_deps}")

        # Find max step number from kept steps
        max_step = 0
        if kept_steps:
            max_step = max(int(step['step']) for step in kept_steps)
        else:
            # If no steps were kept, there might be a problem
            print("WARNING: No steps from original pathway were kept!")

        print(f"Max step number from kept steps: {max_step}")

        # Adjust new steps (renumber them starting from max_step + 1)
        new_steps = []
        step_mapping = {}

        for idx, step in enumerate(new_result.get('steps', [])):

            new_step_num = max_step + 1 + idx
            step_mapping[step['step']] = str(new_step_num)

            adjusted_step = step.copy()
            adjusted_step['step'] = str(new_step_num)
            new_steps.append(adjusted_step)

        print(f"New steps after renumbering: {[s['step'] for s in new_steps]}")

        # Adjust dependencies for new steps
        new_deps = {}
        for old_num, deps in new_result.get('dependencies', {}).items():
            new_num = step_mapping[old_num]
            # Map old step numbers to new step numbers in dependencies
            new_deps[new_num] = [step_mapping[d] for d in deps]

        print(f"New dependencies after renumbering: {new_deps}")

        # Connect the new branch to the left connection if it exists
        if left_connection and new_steps:
            first_new_step = new_steps[0]['step']
            if left_connection in kept_deps:
                kept_deps[left_connection].append(first_new_step)
            else:
                kept_deps[left_connection] = [first_new_step]
            print(
                f"Connected left step {left_connection} to new step {first_new_step}"
            )

        # Merge results - this should include ALL kept steps AND new steps
        merged_steps = kept_steps + new_steps
        merged_deps = {**kept_deps, **new_deps}

        # Verify the merging was successful
        print(f"Merged steps: {[s['step'] for s in merged_steps]}")
        print(f"Merged dependencies: {merged_deps}")

        # Create the final merged result
        merged_result = {'steps': merged_steps, 'dependencies': merged_deps}

        # Final debug check
        if not merged_steps:
            print("ERROR: No steps in final merged result!")

            return jsonify({"error": "No steps in final merged result"}), 500

        # Store the merged result in a JSON file
        save_result(smiles, merged_result)

        print("\n=== Partial Rerun Complete ===")
        return jsonify(merged_result), 200

    except Exception as e:
        print(f"\nERROR in partial rerun:")
        print(f"Exception: {str(e)}")
        print(f"Traceback: {traceback.format_exc()}")
        return jsonify({"error":
                        f"Error in partial retrosynthesis: {str(e)}"}), 500


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)
