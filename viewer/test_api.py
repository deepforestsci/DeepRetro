from flask import Flask, jsonify, request
from flask_cors import CORS
import traceback

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# The JSON data
REACTION_DATA = {
    "dependencies": {
        "1": ["2"],
        "2": ["3"],
        "3": ["4"],
        "4": []
    },
    "steps": [
        {
            "step": "1",
            "reactants": [
                {
                    "smiles": "NC(=O)c1cncc(F)n1",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "C5H4FN3O",
                        "mass": 141.033839968
                    }
                }
            ],
            "reagents": [
                {
                    "smiles": "O",
                    "reagent_metadata": {
                        "name": "",
                        "chemical_formula": "H2O",
                        "mass": 18.010564684
                    }
                }
            ],
            "products": [
                {
                    "smiles": "NC(=O)C1=NC(F)=CN=C1O",
                    "product_metadata": {
                        "name": "",
                        "chemical_formula": "C5H4FN3O2",
                        "mass": 157.028754588
                    }
                }
            ],
            "conditions": [],
            "reactionmetrics": [
                {
                    "scalabilityindex": "10",
                    "confidenceestimate": 0.9,
                    "closestliterature": ""
                }
            ]
        },
        {
            "step": "2",
            "reactants": [
                {
                    "smiles": "O=C(O)O",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "CH2O3",
                        "mass": 62.000393923999994
                    }
                },
                {
                    "smiles": "N#Cc1cncc(F)n1",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "C5H2FN3",
                        "mass": 123.023275284
                    }
                }
            ],
            "reagents": [],
            "products": [
                {
                    "smiles": "NC(=O)c1cncc(F)n1",
                    "product_metadata": {
                        "name": "",
                        "chemical_formula": "C5H4FN3O",
                        "mass": 141.033839968
                    }
                }
            ],
            "conditions": [],
            "reactionmetrics": [
                {
                    "scalabilityindex": "10",
                    "confidenceestimate": 0.99,
                    "closestliterature": ""
                }
            ]
        },
        {
            "step": "3",
            "reactants": [
                {
                    "smiles": "C[Si](C)(C)C#N",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "C4H9NSi",
                        "mass": 99.050425818
                    }
                },
                {
                    "smiles": "[O-][n+]1ccncc1F",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "C4H3FN2O",
                        "mass": 114.022940936
                    }
                }
            ],
            "reagents": [],
            "products": [
                {
                    "smiles": "N#Cc1cncc(F)n1",
                    "product_metadata": {
                        "name": "",
                        "chemical_formula": "C5H2FN3",
                        "mass": 123.023275284
                    }
                }
            ],
            "conditions": [],
            "reactionmetrics": [
                {
                    "scalabilityindex": "8",
                    "confidenceestimate": 0.98,
                    "closestliterature": ""
                }
            ]
        },
        {
            "step": "4",
            "reactants": [
                {
                    "smiles": "Fc1cnccn1",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "C4H3FN2",
                        "mass": 98.028026316
                    }
                },
                {
                    "smiles": "O=C(OO)c1cccc(Cl)c1",
                    "reactant_metadata": {
                        "name": "",
                        "chemical_formula": "C7H5ClO3",
                        "mass": 171.99272169999998
                    }
                }
            ],
            "reagents": [],
            "products": [
                {
                    "smiles": "[O-][n+]1ccncc1F",
                    "product_metadata": {
                        "name": "",
                        "chemical_formula": "C4H3FN2O",
                        "mass": 114.022940936
                    }
                }
            ],
            "conditions": [],
            "reactionmetrics": [
                {
                    "scalabilityindex": "10",
                    "confidenceestimate": 0.74,
                    "closestliterature": ""
                }
            ]
        }
    ]
}

# API key validation (in a production environment, use a more secure method)
API_KEY = "your-secure-api-key"

def validate_api_key():
    api_key = request.headers.get('X-API-KEY')
    return api_key == API_KEY

@app.route('/reactions', methods=['GET', 'POST'])
def handle_reactions():
    try:
        # Validate API key
        if not validate_api_key():
            return jsonify({"error": "Invalid API key"}), 401

        if request.method == 'GET':
            return jsonify(REACTION_DATA)
        
        elif request.method == 'POST':
            data = request.get_json()
            if not data or 'smiles' not in data:
                return jsonify({"error": "Missing SMILES data"}), 400
            
            # For now, return the same reaction data regardless of input SMILES
            # In a real application, you would process the SMILES here
            return jsonify(REACTION_DATA)

    except Exception as e:
        # Log the full error in your production environment
        print(traceback.format_exc())
        return jsonify({"error": str(e)}), 500

@app.errorhandler(404)
def not_found(e):
    return jsonify({"error": "Resource not found"}), 404

@app.errorhandler(405)
def method_not_allowed(e):
    return jsonify({"error": "Method not allowed"}), 405

if __name__ == '__main__':
    app.run(debug=True)