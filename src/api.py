from flask import Flask, request, jsonify
from flask_cors import CORS
from functools import wraps
from dotenv import load_dotenv
from rdkit import Chem

import rootutils

root_dir = rootutils.setup_root(__file__,
                                indicator=".project-root",
                                pythonpath=True)

load_dotenv()
from src.main import main
from src.cache import clear_cache_for_molecule

app = Flask(__name__)
CORS(app,
     resources={
         r"/api/*": {
             "origins":
             "http://ec2-3-142-141-95.us-east-2.compute.amazonaws.com:8000"
         }
     })

# Predefined API key for authentication
API_KEY = "your-secure-api-key"


# Authentication decorator
def require_api_key(f):

    @wraps(f)
    def decorated_function(*args, **kwargs):
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
    result = main(smiles)
    # try:
    #     result = run_prithvi(smiles)
    # except Exception as e:
    #     print(e)
    #     return jsonify({"error": "Error in retrosynthesis"}), 500
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

    molecule = data['molecule']

    # Clear the cache for the molecule
    clear_cache_for_molecule(molecule)

    # Rerun retrosynthesis
    result = main(molecule)
    # try:
    #     result = run_prithvi(smiles)
    # except Exception as e:
    #     print(e)
    #     return jsonify({"error": "Error in retrosynthesis"}), 500
    return jsonify(result), 200


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)
