#!/usr/bin/env python
# Make sure you're running flask.py in a separate terminal first
# Try run  $ ./flask.py "C(C(F)(F)F)O"
# This is only for test. Since the flask.py return images.
# These base64 image data looks like a huge jumble of letters, numbers, and symbols.

import requests
import click
from flask import render_template

url = "http://127.0.0.1:5000/predict"


@click.command()
@click.argument(
    "smiles",
    type=str,
    default="C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O",
)
def predict(smiles):
    """Predict the 19F NMR spectrum of a molecule given its SMILES string"""
    headers = {"Accept": "application/json"}  # Explicitly request JSON
    response = requests.get(f"{url}/{smiles}", headers=headers)
    # Check if the response content is valid

    if response.status_code == 200:
        try:
            result = (
                response.json()
            )  # This may throw an error if the response isn't valid JSON
            plot_data = result.get("plot_data", "No plot data available")
            table_data = result.get("table_data", "No table data available")
            structure_image_base64 = result.get(
                "structure_image_base64", "No image available"
            )
            # Print results for CLI usage
            click.echo("Prediction Results:")
            click.echo(f"Plot Data: {plot_data}")
            click.echo(f"Table Data: {table_data}")
            click.echo(f"Structure Image: {structure_image_base64}")
        except ValueError as e:
            click.echo(f"Error: Response is not valid JSON - {e}")
    else:
        click.echo(f"Error: {response.status_code} - {response.text}")


if __name__ == "__main__":
    predict()
