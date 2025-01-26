from flask import Flask, request, render_template
from predictor import predictor  # Make sure to import the predictor function

app = Flask(__name__)  # Corrected: use __name__ instead of name


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        smiles = request.form.get("smiles")  # Get the SMILES input from the form

        if smiles:
            # Call the predictor function with the SMILES input
            # results = predictor(smiles)
            plot_data, table_data, structure_image_base64 = predictor(smiles)
            return render_template(
                "index.html",
                plot_data=plot_data,
                table_data=table_data,
                structure_image_base64=structure_image_base64,
            )
        else:
            return "Please provide a valid SMILES string."
    return render_template("index.html")


if __name__ == "__main__":  # Corrected: use __name__ == "__main__"
    app.run(debug=True, host="0.0.0.0", port=5002)
