from flask import Flask, request, jsonify
from flask_cors import CORS
from functools import wraps
from dotenv import load_dotenv
from rdkit import Chem

import rootutils
import os
import json
import traceback

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

load_dotenv()
from src.main import main
from src.cache import clear_cache_for_molecule
from src.variables import AZ_MODEL_LIST

# Load advanced settings config once at startup
config_path = os.path.join(root_dir, 'config', 'advanced_settings.json')
with open(config_path) as f:
    advanced_config = json.load(f)

app = Flask(__name__)
CORS(app)

# API key loaded from environment variable
API_KEY = os.getenv('API_KEY')

# For testing: Uncomment the line below to see the loaded API_KEY
# print(f"DEBUG: API_KEY loaded from environment: {API_KEY}")

if not API_KEY:
    print("CRITICAL ERROR: The 'API_KEY' environment variable is not set.")
    print(
        "Please set this variable in your .env file or your system environment."
    )
    print("Example: API_KEY='your-chosen-secret-key'")
    exit(1)  # Exit if the API key is not configured

# File path for storing results
PARTIAL_JSON_PATH = "partial.json"


# Functions for JSON file storage
def save_result(smiles, result):
    """Save retrosynthesis result to partial.json file."""
    # Create a data structure with the SMILES and result
    data = {"smiles": smiles, "result": result}

    # Save the result to the file (overwriting any existing data)
    with open(PARTIAL_JSON_PATH, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"Result for {smiles} saved to {PARTIAL_JSON_PATH}")
    return PARTIAL_JSON_PATH


def load_result():
    """Load the most recent retrosynthesis result from partial.json."""
    # Check if the file exists
    if not os.path.exists(PARTIAL_JSON_PATH):
        print(f"No result file found at {PARTIAL_JSON_PATH}")
        return None, None

    # Load the result from the file
    try:
        with open(PARTIAL_JSON_PATH, 'r') as f:
            data = json.load(f)

        smiles = data.get("smiles")
        result = data.get("result")

        print(f"Result for {smiles} loaded from {PARTIAL_JSON_PATH}")
        return smiles, result
    except Exception as e:
        print(f"Error loading result file: {e}")
        return None, None


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
            "error": "SMILES string is required. Please include a 'smiles' field"
        }), 400

    smiles = data['smiles']
    if not Chem.MolFromSmiles(smiles):
        return jsonify({"error": "Invalid SMILES string"}), 400

    defaults = advanced_config['defaults']
    model_type = data.get('model_type', defaults['model_type'])
    if model_type not in advanced_config['llm_models']:
        return jsonify({"error": f"Unsupported model type: {model_type}"}), 400
    model_config = advanced_config['llm_models'][model_type]
    llm = model_config['internal_name']

    advanced_prompt = data.get('advanced_prompt', defaults['advanced_prompt'])
    if isinstance(advanced_prompt, str):
        advanced_prompt = advanced_prompt.lower() == "true"
    if advanced_prompt and not model_config['supports_advanced_prompt']:
        advanced_prompt = False
    if advanced_prompt:
        llm += ":adv"

    az_model = data.get('model_version', defaults['model_version'])
    if az_model not in advanced_config['az_models']:
        return jsonify({"error": f"Unsupported AZ model: {az_model}"}), 400

    stability_flag = data.get('stability_flag', defaults['stability_flag'])
    if isinstance(stability_flag, str):
        stability_flag = stability_flag.lower() == "true"
    if stability_flag and not model_config['supports_stability_check']:
        stability_flag = False

    hallucination_check = data.get('hallucination_check', defaults['hallucination_check'])
    if isinstance(hallucination_check, str):
        hallucination_check = hallucination_check.lower() == "true"
    if hallucination_check and not model_config['supports_hallucination_check']:
        hallucination_check = False

    try:
        result = main(
            smiles=smiles,
            llm=llm,
            az_model=az_model,
            stability_flag=str(stability_flag),
            hallucination_check=str(hallucination_check)
        )
        save_result(smiles, result)
    except Exception as e:
        print(e)
        return jsonify({"error": f"Error in retrosynthesis: {str(e)}. Please rerun."}), 500
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

    # Check if the SMILES string is valid
    if not Chem.MolFromSmiles(molecule):
        return jsonify({"error": "Invalid SMILES string"}), 400

    # -----------------
    # Advanced model - DeepSeek-R1
    model_type = "claude3"  # Default is Claude 3 Opus
    try:
        if 'model_type' in data:
            model_type = data['model_type'].lower()
            assert model_type in ["claude3", "claude37", "deepseek"]
            print(f"USING MODEL TYPE: {model_type}")
    except Exception as e:
        print(f"Error processing model type: {e}")
        model_type = "claude3"
        print(f"FALLING BACK TO DEFAULT MODEL TYPE: {model_type}")

    # Select the appropriate LLM based on model_type
    if model_type == "deepseek":
        llm = "fireworks_ai/accounts/fireworks/models/deepseek-r1"
    elif model_type == "claude37":
        llm = "anthropic/claude-3-7-sonnet-20250219"
    elif model_type == "claude4opus":
        llm = "anthropic/claude-opus-4-20250514"
    elif model_type == "claude4sonnet":
        llm = "anthropic/claude-sonnet-4-20250514"
    else:  # Default to Claude 3 Opus
        llm = "claude-3-opus-20240229"

    # -----------------
    # Advanced prompt handling
    advanced_prompt = False
    try:
        if 'advanced_prompt' in data:
            advanced_prompt = data['advanced_prompt']
            if advanced_prompt.lower() == "true":
                advanced_prompt = True
    except Exception as e:
        print(f"Error processing advanced prompt: {e}")
        advanced_prompt = False

    if advanced_prompt:
        llm = llm + ":adv"

    # -----------------
    # Choose AiZynthFinder model
    az_model = "USPTO"
    try:
        if 'model_version' in data:
            az_model = data['model_version']
            assert az_model in AZ_MODEL_LIST
    except Exception as e:
        print(f"Error processing model version: {e}")
        az_model = "USPTO"

    # -----------------
    # Stability check flag
    stability_flag = "False"
    try:
        if 'stability_flag' in data:
            stability_flag = data['stability_flag']
            assert stability_flag.lower() in ["false", "true"]
    except Exception as e:
        print(f"Error processing stability flag: {e}")
        stability_flag = "False"

    # -----------------
    # Hallucination check flag
    hallucination_check = "False"
    try:
        if 'hallucination_check' in data:
            hallucination_check = data['hallucination_check']
            assert hallucination_check.lower() in ["false", "true"]
    except Exception as e:
        print(f"Error processing hallucination check: {e}")
        hallucination_check = "False"

    # -----------------
    # Rerun retrosynthesis
    try:
        result = main(smiles=molecule,
                      llm=llm,
                      az_model=az_model,
                      stability_flag=stability_flag,
                      hallucination_check=hallucination_check)

        # Store the result in partial.json
        save_result(molecule, result)
    except Exception as e:
        print(e)
        return jsonify(
            {"error":
             f"Error in retrosynthesis: {str(e)}. Please rerun."}), 500

    return jsonify(result), 200


@app.route('/api/partial_rerun', methods=['POST'])
@require_api_key
def partial_rerun():
    """
    Endpoint to partially rerun retrosynthesis from a specific step.
    Uses the stored results from the most recent retrosynthesis run in partial.json.
    
    When rerunning a step, we remove that step and everything to its right in the synthesis pathway.
    """
    print(
        "\n===================== STARTING PARTIAL RERUN PROCESS ====================="
    )

    data = request.get_json()
    print(f"RECEIVED REQUEST DATA: {json.dumps(data, indent=2)}")

    try:
        smiles = data['smiles']
        print(f"TARGET MOLECULE: {smiles}")
        from_step = int(data['steps'])
        print(f"TARGET STEP FOR RERUN: {from_step}")

        # Load previous results from partial.json
        stored_smiles, original_result = load_result()
        print(f"LOADED FROM STORAGE - SMILES: {stored_smiles}")

        if not original_result:
            print("ERROR: No results found in storage")
            return jsonify({
                "error":
                "No previous results found. Run retrosynthesis first."
            }), 400

        if stored_smiles != smiles:
            print(
                f"ERROR: Stored SMILES ({stored_smiles}) doesn't match requested SMILES ({smiles})"
            )
            return jsonify({
                "error":
                "No results found for this molecule. Run retrosynthesis first."
            }), 400

        print(f"FOUND STORED RESULT FOR SMILES: {smiles}")

        # Print the original steps for debugging
        print(f"ORIGINAL STEPS:")
        for step in original_result.get('steps', []):
            print(f"  Step {step.get('step', 'unknown')}: {json.dumps(step)}")

        # Print the original dependency structure for debugging
        print(
            f"ORIGINAL DEPENDENCIES: {json.dumps(original_result.get('dependencies', {}), indent=2)}"
        )

        # Get the step we're rerunning from the original result
        target_step = next((step for step in original_result['steps']
                            if int(step['step']) == from_step), None)

        if not target_step:
            print(f"ERROR: Step {from_step} not found in synthesis pathway")
            return jsonify({
                "error":
                f"Step {from_step} not found in the synthesis pathway"
            }), 404

        # Print the target step for debugging
        print(f"TARGET STEP DETAILS: {json.dumps(target_step, indent=2)}")

        # Get the target molecule from the step's product
        # This is what we want to resynthesize with a new approach
        if not target_step.get('products') or not isinstance(
                target_step['products'], list) or len(
                    target_step['products']) == 0:
            print(f"ERROR: Step {from_step} doesn't have valid products")
            return jsonify(
                {"error":
                 f"Step {from_step} doesn't have valid products"}), 400

        # Log all products in the target step
        print(f"PRODUCTS IN TARGET STEP:")
        for i, product in enumerate(target_step['products']):
            print(f"  Product {i}: {json.dumps(product)}")

        # Get the first product's SMILES (this is what we'll resynthesize)
        start_molecule = target_step['products'][0]['smiles']
        print(f"\nSTARTING NEW SYNTHESIS FROM MOLECULE: {start_molecule}")

        # -----------------
        # Advanced model - DeepSeek-R1
        model_type = "claude3"  # Default is Claude 3 Opus
        try:
            if 'model_type' in data:
                model_type = data['model_type'].lower()
                assert model_type in ["claude3", "claude37", "deepseek"]
                print(f"USING MODEL TYPE: {model_type}")
        except Exception as e:
            print(f"Error processing model type: {e}")
            model_type = "claude3"
            print(f"FALLING BACK TO DEFAULT MODEL TYPE: {model_type}")

        # Select the appropriate LLM based on model_type
        if model_type == "deepseek":
            llm = "fireworks_ai/accounts/fireworks/models/deepseek-r1"
        elif model_type == "claude37":
            llm = "anthropic/claude-3-7-sonnet-20250219"
        elif model_type == "claude4opus":
            llm = "anthropic/claude-opus-4-20250514"
        elif model_type == "claude4sonnet":
            llm = "anthropic/claude-sonnet-4-20250514"
        else:  # Default to Claude 3 Opus
            llm = "claude-3-opus-20240229"
        print(f"SELECTED LLM: {llm}")

        # -----------------
        # Advanced prompt handling
        advanced_prompt = False
        try:
            if 'advanced_prompt' in data:
                advanced_prompt = data['advanced_prompt']
                if advanced_prompt.lower() == "true":
                    advanced_prompt = True
                print(f"USING ADVANCED PROMPT: {advanced_prompt}")
        except Exception as e:
            print(f"ERROR PROCESSING ADVANCED PROMPT: {e}")
            advanced_prompt = False

        if advanced_prompt:
            llm = llm + ":adv"
            print(f"UPDATED LLM WITH ADVANCED PROMPT: {llm}")

        # -----------------
        # Choose AiZynthFinder model
        az_model = "USPTO"
        try:
            if 'model_version' in data:
                az_model = data['model_version']
                assert az_model in AZ_MODEL_LIST
                print(f"USING MODEL VERSION: {az_model}")
        except Exception as e:
            print(f"ERROR PROCESSING MODEL VERSION: {e}")
            az_model = "USPTO"
            print(f"FALLING BACK TO DEFAULT MODEL: {az_model}")

        # -----------------
        # Stability check flag
        stability_flag = "False"
        try:
            if 'stability_flag' in data:
                stability_flag = data['stability_flag']
                assert stability_flag.lower() in ["false", "true"]
                print(f"USING STABILITY FLAG: {stability_flag}")
        except Exception as e:
            print(f"ERROR PROCESSING STABILITY FLAG: {e}")
            stability_flag = "False"
            print(f"FALLING BACK TO DEFAULT STABILITY FLAG: {stability_flag}")

        # -----------------
        # Hallucination check flag
        hallucination_check = "False"
        try:
            if 'hallucination_check' in data:
                hallucination_check = data['hallucination_check']
                assert hallucination_check.lower() in ["false", "true"]
                print(f"USING HALLUCINATION CHECK: {hallucination_check}")
        except Exception as e:
            print(f"ERROR PROCESSING HALLUCINATION CHECK: {e}")
            hallucination_check = "False"
            print(
                f"FALLING BACK TO DEFAULT HALLUCINATION CHECK: {hallucination_check}"
            )

        # Run new synthesis on the starting molecule
        print(f"\nCALLING MAIN FUNCTION WITH PARAMETERS:")
        print(f"  SMILES: {start_molecule}")
        print(f"  LLM: {llm}")
        print(f"  AZ MODEL: {az_model}")
        print(f"  STABILITY FLAG: {stability_flag}")
        print(f"  HALLUCINATION CHECK: {hallucination_check}")

        try:
            new_result = main(smiles=start_molecule,
                              llm=llm,
                              az_model=az_model,
                              stability_flag=stability_flag,
                              hallucination_check=hallucination_check)
            print(
                f"NEW RETROSYNTHESIS RESULT: {json.dumps(new_result, indent=2)}"
            )

            # Note if the retrosynthesis result is empty, but treat it as a valid result
            if not new_result.get('steps'):
                print(
                    f"NOTE: Empty synthesis result for {start_molecule}. This molecule doesn't require further decomposition or is a building block."
                )

                # Initialize with empty steps and dependencies if not present
                if 'steps' not in new_result:
                    new_result['steps'] = []
                    print("INITIALIZED empty steps list")
                if 'dependencies' not in new_result:
                    new_result['dependencies'] = {}
                    print("INITIALIZED empty dependencies dict")

        except Exception as e:
            print(
                f"ERROR RUNNING RETROSYNTHESIS ON MOLECULE {start_molecule}:")
            print(f"  EXCEPTION: {str(e)}")
            print(f"  TRACEBACK: {traceback.format_exc()}")
            return jsonify({
                "error":
                f"Error running retrosynthesis on {start_molecule}: {str(e)}"
            }), 500

        # The steps to remove are the target step and everything it depends on
        # In other words, the target step and everything to its right in the synthesis pathway
        steps_to_remove = {str(from_step)}
        steps_to_check = [str(from_step)]

        print(f"\nIDENTIFYING STEPS TO REMOVE:")
        print(f"  STARTING WITH STEP: {from_step}")

        # Identify all steps that the target step depends on (to the right in the pathway)
        # This is the set of steps that will be replaced by the new synthesis
        while steps_to_check:
            current_step = steps_to_check.pop(0)
            print(f"  CHECKING DEPENDENCIES FOR STEP: {current_step}")
            if current_step in original_result.get('dependencies', {}):
                print(
                    f"    FOUND DEPENDENCIES: {original_result['dependencies'][current_step]}"
                )
                for dep_step in original_result['dependencies'][current_step]:
                    if dep_step not in steps_to_remove:
                        steps_to_remove.add(dep_step)
                        steps_to_check.append(dep_step)
                        print(f"    ADDED STEP {dep_step} TO REMOVAL LIST")

        print(
            f"STEPS TO REMOVE (TARGET AND EVERYTHING TO ITS RIGHT): {steps_to_remove}"
        )

        # We need to identify what step the target step is connected to on its left
        # This is where we'll connect the new synthesis
        left_connection = None
        for step, deps in original_result.get('dependencies', {}).items():
            if str(from_step) in deps and step not in steps_to_remove:
                left_connection = step
                print(
                    f"FOUND LEFT CONNECTION: STEP {step} DEPENDS ON STEP {from_step}"
                )
                break

        print(
            f"LEFT CONNECTION (STEP THAT THE TARGET CONNECTS TO ON THE LEFT): {left_connection}"
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

        print(f"KEPT STEPS: {[s['step'] for s in kept_steps]}")
        print(f"KEPT DEPENDENCIES: {kept_deps}")

        # Find max step number from kept steps
        max_step = 0
        if kept_steps:
            max_step = max(int(step['step']) for step in kept_steps)
            print(f"MAX STEP NUMBER FROM KEPT STEPS: {max_step}")
        else:
            # If no steps were kept, there might be a problem
            print("WARNING: No steps from original pathway were kept!")

        # Adjust new steps (renumber them starting from max_step + 1)
        new_steps = []
        step_mapping = {}

        print(f"\nRENUMBERING NEW STEPS STARTING FROM {max_step + 1}:")
        for idx, step in enumerate(new_result.get('steps', [])):
            new_step_num = max_step + 1 + idx
            old_step_num = str(
                step['step'])  # Convert to string for consistency
            step_mapping[old_step_num] = str(new_step_num)
            print(f"  MAPPING STEP {old_step_num} -> {new_step_num}")

            adjusted_step = step.copy()
            adjusted_step['step'] = str(new_step_num)
            new_steps.append(adjusted_step)

        print(f"NEW STEPS AFTER RENUMBERING: {[s['step'] for s in new_steps]}")

        # Adjust dependencies for new steps
        new_deps = {}
        print(f"\nADJUSTING DEPENDENCIES FOR NEW STEPS:")
        for old_num, deps in new_result.get('dependencies', {}).items():
            old_num = str(old_num)  # Convert to string for consistency
            new_num = step_mapping[old_num]
            print(f"  STEP {old_num} -> {new_num} DEPENDS ON:")

            # Map old step numbers to new step numbers in dependencies
            mapped_deps = []
            for d in deps:
                d = str(d)  # Convert to string for consistency
                mapped_d = step_mapping[d]
                mapped_deps.append(mapped_d)
                print(f"    OLD DEP {d} -> NEW DEP {mapped_d}")

            new_deps[new_num] = mapped_deps

        print(f"NEW DEPENDENCIES AFTER RENUMBERING: {new_deps}")

        # Connect the new branch to the left connection if it exists
        if left_connection and new_steps:
            first_new_step = new_steps[0]['step']
            print(
                f"\nCONNECTING LEFT STEP {left_connection} TO NEW STEP {first_new_step}"
            )

            if left_connection in kept_deps:
                kept_deps[left_connection].append(first_new_step)
                print(
                    f"  ADDED {first_new_step} TO EXISTING DEPENDENCIES OF {left_connection}: {kept_deps[left_connection]}"
                )
            else:
                kept_deps[left_connection] = [first_new_step]
                print(
                    f"  CREATED NEW DEPENDENCIES FOR {left_connection}: {kept_deps[left_connection]}"
                )

        # Merge results - this should include ALL kept steps AND new steps
        merged_steps = kept_steps + new_steps
        merged_deps = {**kept_deps, **new_deps}

        # Verify the merging was successful
        print(f"\nFINAL MERGED STEPS: {[s['step'] for s in merged_steps]}")
        print(f"FINAL MERGED DEPENDENCIES: {merged_deps}")

        # Final debug check - NOTE: Empty results are now considered valid when a molecule
        # is identified as a building block or starting material
        if not merged_steps and not new_result.get('steps'):
            print(
                "NOTE: Final result has no steps. The molecule was identified as a building block or starting material."
            )
            # Create a minimal result that acknowledges this molecule as a building block
            merged_result = {
                'smiles':
                smiles,  # Add the original SMILES
                'steps': [],
                'dependencies': {},
                'message':
                f"The molecule {start_molecule} was identified as a building block or starting material and doesn't require further synthesis."
            }
        elif not merged_steps:
            print("ERROR: NO STEPS IN FINAL MERGED RESULT!")
            return jsonify({"error": "No steps in final merged result"}), 500
        else:
            # Create the normal merged result
            merged_result = {
                'smiles': smiles,  # Add the original SMILES to the result
                'steps': merged_steps,
                'dependencies': merged_deps
            }

        # Store the merged result in partial.json
        save_result(smiles, merged_result)
        print(f"SAVED MERGED RESULT TO {PARTIAL_JSON_PATH}")

        print(
            "\n===================== PARTIAL RERUN COMPLETE ====================="
        )
        return jsonify(merged_result), 200

    except Exception as e:
        print(f"\nFATAL ERROR IN PARTIAL RERUN:")
        print(f"  EXCEPTION: {str(e)}")
        print(f"  TRACEBACK: {traceback.format_exc()}")
        return jsonify({"error":
                        f"Error in partial retrosynthesis: {str(e)}"}), 500


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)
