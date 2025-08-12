#!/usr/bin/env python
# In terminal, run $ ./flask.py
# Open a web browser and go to http://127.0.0.1:5000
# Follow the instructions to make a prediction.

from flask import Flask, request, jsonify, render_template
from predictor import predictor

app = Flask(__name__)

@app.route('/')
def home():
    return 'API is running! Add /predict/SMILES to the URL to make a prediction.  ' \
           'Example: /predict/C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O'

@app.route('/predict', methods=['GET'])
def index():
    return render_template("index.html")

@app.route('/predict/<string:smiles>', methods=['GET', 'POST'])
def predict(smiles):
    """Predict the 19F NMR spectrum of a molecule given its SMILES string"""
    if not smiles:
            return jsonify({"error": "Please provide a SMILES string."})
    
    try:
        # Call the predictor function with the SMILES input
        plot_data, table_data, structure_image_base64 = predictor(smiles)

        # When a web browser (or any client expecting HTML) makes a GET request, it renders index.html.
        if "text/html" in request.accept_mimetypes:
            return render_template(
                "index_rest.html",
                plot_data=plot_data,
                table_data=table_data,
                structure_image_base64=structure_image_base64,
            )
        
        # When a CLI or API client (like requests.get()) makes a GET request, it returns JSON.
        elif "application/json" in request.headers.get("Accept", ""):
            return jsonify({
                "plot_data": plot_data,
                "table_data": table_data,
                "structure_image_base64": structure_image_base64
            }), 200
        # return response, status code, headers

    except ValueError as e:
        # Pass error message to the template
        return jsonify({"error": str(e)}), 400
    
    except Exception as e:
            # Handle unexpected exceptions
        return jsonify({"error": f"Unexpected error: {str(e)}"}), 500

if __name__ == "__main__":
    app.run(debug = True)