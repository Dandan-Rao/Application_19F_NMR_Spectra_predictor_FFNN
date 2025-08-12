#!/usr/bin/env python
from flask import Flask, request, render_template
from predictor import predictor  # Make sure to import the predictor function

app = Flask(__name__)  # Corrected: use __name__ instead of name


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        smiles = request.form.get("smiles")  # Get the SMILES input from the form

        if not smiles:
            return render_template(
                "index.html", error_message="Please provide a valid SMILES string."
            )

        try:
            # Call the predictor function with the SMILES input
            plot_data, table_data, structure_image_base64 = predictor(smiles)
            return render_template(
                "index.html",
                plot_data=plot_data,
                table_data=table_data,
                structure_image_base64=structure_image_base64,
            )
        except ValueError as e:
            # Pass error message to the template
            return render_template("index.html", error_message=str(e))
        except Exception as e:
            # Handle unexpected exceptions
            return render_template(
                "index.html", error_message=f"Unexpected error: {str(e)}"
            )

    return render_template("index.html")


if __name__ == "__main__":  # Corrected: use __name__ == "__main__"
    app.run(debug=False, host="0.0.0.0", port=5002)
