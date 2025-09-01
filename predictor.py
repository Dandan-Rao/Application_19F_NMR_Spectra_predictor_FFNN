# +
import os
import io
import sys
import base64
import pandas as pd

import matplotlib

matplotlib.use("Agg")  # Set non-GUI backend before importing pyplot
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw
from utils import *

from tensorflow.keras.models import load_model

pd.options.mode.chained_assignment = None  # Suppress the SettingWithCopyWarning
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def get_sdf_file(smiles):
    """
    Create .sdf file and temporarily save as /temp/temp.sdf
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid SMILES string")

    if not any(atom.GetSymbol() == "F" for atom in mol.GetAtoms()):
        raise ValueError("No F in the molecule")

    # Transform SMILES to canonical SMILES
    smiles = Chem.MolToSmiles(mol)

    mol = Chem.MolFromSmiles(smiles)
    sdf = Transform_SMILE_to_3D_conformation(smiles)

    # Save .sdf file
    filename = os.path.join("temp", "temp.sdf")

    w = Chem.SDWriter(filename)
    w.write(sdf)
    w.close()

    # Create pictures for the molecule
    for _, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))

        # Save the molecule image with the Code as part of the file name
        file_path = os.path.join("temp", "temp.png")
        Draw.MolToFile(mol, file_path, size=(600, 400))

# +
def get_test_fluorianted_compounds_info(smiles, train_dataset):
    """
    Create a dataframe containing Code, SMILES, and atom indexs as columns.
    """
    if smiles in train_dataset["SMILES"].values:
        dataset = train_dataset[train_dataset["SMILES"] == smiles]
        if len(dataset) > 1:
            dataset = dataset.iloc[[0], :]
    #         dataset['Code'] = 'temp'

    else:
        # Define columns from '0' to '70' and add an additional column 'Code'
        columns = ["Code", "SMILES"] + [str(i) for i in range(71)]

        # Create the DataFrame with 71 columns filled with None and 'Code' filled with 'temp'
        dataset = pd.DataFrame([["temp"] + [smiles] + [None] * 71], columns=columns)
    return dataset



def get_HOSE_prediction_results_table(
    HOSE_Code_database_file_path, test_fluorinated_compounds
):
    HOSE_Code_database = pd.read_csv(HOSE_Code_database_file_path)
    # Transform column names to int where possible
    HOSE_Code_database.columns = [
        convert_column_name(name) for name in HOSE_Code_database.columns
    ]

    HOSE_codes_test = getHoseCodeContent(test_fluorinated_compounds)

    # Get HOSE Code and corresponding 19F NMR values using train dataset
    sphere_dics = getTrainDictionary_HOSE(HOSE_Code_database)

    HOSE_Code_database["NMR_Peaks"] = HOSE_Code_database["NMR_Peaks"].apply(
        pd.to_numeric, errors="coerce"
    )

    # Get prediction results and corresponding similarity levels for the validation dataset
    prediction, similarity_levels = HOSE_Model(
        sphere_dics, HOSE_codes_test, HOSE_Code_database["NMR_Peaks"].mean()
    )
    # Validation dataset
    results = getResults_HOSE(prediction, similarity_levels, HOSE_codes_test)

    return results


def safe_split(index_value):
    parts = index_value.split("_", 1)
    if len(parts) == 2:
        return parts
    elif len(parts) == 1:
        return [parts[0], None]  # Only one part, return None for the second
    else:
        return [None, None]  # No parts, return None for both


def display_results(results):
    try:
        # Set style for better visualization
        plt.style.use("seaborn-v0_8-darkgrid")

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), width_ratios=[2, 1])
        fig.subplots_adjust(hspace=0.1)

        # Prepare data
        real_PFAS_spectra_df = results[
            ["fluorinated_compounds", "atom_index", "actual"]
        ]
        code = real_PFAS_spectra_df["fluorinated_compounds"].iloc[0]
        real_PFAS_spectra_df = real_PFAS_spectra_df.pivot(
            index="fluorinated_compounds", columns="atom_index", values="actual"
        )

        predicted_PFAS_spectra_df = results[
            ["fluorinated_compounds", "atom_index", "ensembled_model"]
        ]
        predicted_PFAS_spectra_df = predicted_PFAS_spectra_df.pivot(
            index="fluorinated_compounds",
            columns="atom_index",
            values="ensembled_model",
        )

        # Colors for better contrast
        actual_color = "g"  # Green
        predict_color = "C9"  # Blue

        # Get actual and predicted value counts
        actual = real_PFAS_spectra_df.loc[code, :].value_counts()
        predict = predicted_PFAS_spectra_df.loc[code, :].value_counts()

        # Plot in the first subplot (ax1)
        if actual.empty:
            plt.vlines(
                predict.index, ymin=0, ymax=-predict.values, color="w", label="Actual"
            )

        else:
            ax1.vlines(
                actual.index,
                ymin=0,
                ymax=-actual.values,
                color=actual_color,
                label="Actual",
                linewidth=2,
            )

        # Plot predicted values
        ax1.vlines(
            predict.index,
            ymin=0,
            ymax=predict.values,
            color=predict_color,
            label="Prediction",
            linewidth=2,
        )

        # Add similarity levels
        #     temp = results[results['fluorinated_compounds'] == 'temp']
        for i, j in zip(predict.index, predict.values):
            similarity_level = results[results["ensembled_model"] == i][
                "similarity_levels"
            ]
            if not similarity_level.empty:
                ax1.text(
                    i,
                    j + 0.05,
                    similarity_level.iloc[0],
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="black",
                    bbox=dict(facecolor="white", alpha=0.3),
                )

        # Set plot limits and labels
        if not actual.empty:
            x_min = min(actual.index.min(), predict.index.min()) - 10
            x_max = max(actual.index.max(), predict.index.max()) + 10
            y_max = max(actual.values.max(), predict.values.max()) + 0.2
        else:
            x_min = predict.index.min() - 10
            x_max = predict.index.max() + 10
            y_max = predict.values.max() + 0.2
            ax1.text(
                (x_min + x_max) / 2,
                -y_max / 2,
                "No report value",
                color=actual_color,
                ha="center",
                va="center",
                fontsize=17,
            )

        ax1.set_xlim([x_min, x_max])
        ax1.set_ylim([-y_max, y_max])

        # Add labels
        ax1.text(
            x_min - 2,
            y_max / 2,
            "Prediction",
            color=predict_color,
            rotation=90,
            ha="center",
            va="center",
            fontsize=14,
        )
        ax1.text(
            x_min - 2,
            -y_max / 2,
            "Report",
            color=actual_color,
            rotation=90,
            ha="center",
            va="center",
            fontsize=14,
        )

        ax1.set_xlabel(r"$^{19}$F NMR shift", fontsize=12)
        ax1.set_yticks([])
        ax1.axhline(0, color="k", linewidth=0.5)
        ax1.set_title("NMR Shift Prediction Results", fontsize=14, pad=20)

        # Add confidence level information in the second subplot (ax2)
        confidence_data = {
            "Level": [6, 5, 4, 3, 2, 1],
            "Error": [1.0, 1.8, 4.9, 7.6, 13.7, 13.8],
        }

        ax2.axis("off")
        ax2.set_xlim(0, 2)  # Adjust limits as needed
        ax2.set_ylim(0, 2)
        table_text = (
            "Note:\n"
            + "The number above prediction results indicates\n"
            + "the confidence level of the prediction.\n\n"
            + "Confidence Levels (75% prediction error):\n\n"
        )
        for level, error in zip(confidence_data["Level"], confidence_data["Error"]):
            table_text += f"Level {level}: ±{error} ppm\n"
        ax2.text(
            0,
            1,
            table_text,
            ha="left",
            va="center",
            fontsize=12,
            bbox=dict(facecolor="white", edgecolor="white", alpha=0.9),
        )

        plot_img = io.BytesIO()
        plt.savefig(plot_img, format="png")
        plt.close(fig)  # Close figure to free memory
        plot_img.seek(0)
        plot_base64 = base64.b64encode(plot_img.getvalue()).decode("utf-8")

        details = results[
            [
                "atom_index",
                "similarity_levels",
                "actual",
                "ensembled_model",
                "ensembled_model_error",
            ]
        ].rename(
            columns={
                "atom_index": "Atom Index",
                "similarity_levels": "Confidence Level",
                "actual": "Report Values",
                "ensembled_model": "Prediction Results",
                "ensembled_model_error": "Prediction Error",
            }
        )

        details.reset_index(drop=True, inplace=True)

        # Save the table for web use
        table_data = details.to_html(classes="table table-striped", index=False)

        # Check if the molecular structure image exists and encode it if necessary
        structure_image_path = os.path.join("temp", "temp.png")
        if os.path.exists(structure_image_path):
            with open(structure_image_path, "rb") as structure_file:
                structure_base64 = base64.b64encode(structure_file.read()).decode(
                    "utf-8"
                )
        else:
            structure_base64 = None

    except Exception as e:
        print(f"Error generating structure image: {str(e)}")

    return plot_base64, table_data, structure_base64


# -


def predictor(
    smiles,
    # experimental_data_file_name, 
    train_fluorinated_compounds_file_path=os.path.join("artifacts", "Processed_PFAS_19F_NMR_spectra_data.csv"),
    FFNN_2D_model_path=os.path.join("artifacts", "Application_2D_FFNN_sphere5.h5"), 
    scaler_path=os.path.join("artifacts", "Application_2D_FFNN_scaler_2d_sphere5.pkl"),
    imputer_path=os.path.join("artifacts", "Application_2D_FFNN_imputer_2d_sphere5.pkl"), 
    columns_path=os.path.join("artifacts", "Application_2D_FFNN_column_names_2d_sphere5.pkl"),
    HOSE_Code_database_file_path=os.path.join("artifacts", "HOSE_database_all_fluorianted_compounds.csv"),
):
    try:
            # Load preprocessing tools
        with open(scaler_path, "rb") as file:
            scaler = pickle.load(file)
        with open(imputer_path, "rb") as file:
            imputer = pickle.load(file)
        with open(columns_path, "rb") as file:
            train_columns = pickle.load(file)
        best_model = load_model(FFNN_2D_model_path, safe_mode=False)
        train_dataset = pd.read_csv(train_fluorinated_compounds_file_path, index_col=0)
        # Transform ionic SMILES to neutral SMILES, and transform SMILES to canonical SMILES
        smiles = ionic_to_neutral_smiles(smiles)

        # Generate sdf file from SMILES
        get_sdf_file(smiles)

        dataset = get_test_fluorianted_compounds_info(smiles, train_dataset)

        # # Generate CDK descriptors and Neighbors information
        # get_descriptors_and_neighbors_info()

        # Generate a features table. If SMILES in train dataset then return rows from train, else return empty DataFrame
        dataset = get_test_fluorianted_compounds_info(smiles, train_dataset)

        # Fill the dataset with 2D descriptors if it's empty
        get_2d_descriptors = getAtomicDescriptorsFrom2DNeighbors()
        content = get_2d_descriptors.getDescriptorsFromDataset(dataset=dataset, num_spheres=5)

        content = content.apply(safe_to_numeric)
        content.columns = content.columns.astype(str)

        # NMR_peaks_with_desc = get_features_table(dataset)

        # Get Prediction results from HOSE model
        HOSE_results = get_HOSE_prediction_results_table(
            HOSE_Code_database_file_path, dataset
        )

        y = content["NMR_Peaks"]
        X = content.drop(["NMR_Peaks"], axis=1)
        X = X[train_columns]  # Ensure correct column alignment

         # Impute and scale
        X_imputed = imputer.transform(X)
        X_imputed_df = pd.DataFrame(X_imputed, columns=train_columns, index=X.index)

        X_scaled = scaler.transform(X_imputed_df)
        X_scaled_df = pd.DataFrame(X_scaled, columns=train_columns, index=X.index)

        # FFNN prediction
        FFNN_results_table = get_results_table(best_model=best_model, X=X_scaled_df, y=y)

        # HOSE prediction
        HOSE_results = get_HOSE_prediction_results_table(HOSE_Code_database_file_path, dataset)

        # Combine results
        combined_prediction = HOSE_results.copy()
        combined_prediction.rename(columns={"prediction": "HOSE_model_prediction"}, inplace=True)
        combined_prediction = combined_prediction[["actual", "similarity_levels", "HOSE_model_prediction"]]
        combined_prediction["FFNN_model_prediction"] = FFNN_results_table["prediction"]
        # combined_prediction["FFNN_model_error"] = np.abs(combined_prediction["FFNN_model_prediction"] - combined_prediction["actual"])

        # Ensemble logic
        ensembled_predictions = []
        for _, row in combined_prediction.iterrows():
            if row["similarity_levels"] > 4:
                ensembled_predictions.append(row["HOSE_model_prediction"])
            else:
                ensembled_predictions.append(row["FFNN_model_prediction"])
        combined_prediction["ensembled_model"] = ensembled_predictions

        # Add atom index and compound name
        split_values = [safe_split(idx) for idx in combined_prediction.index]
        combined_prediction["atom_index"] = [val[0] for val in split_values]
        combined_prediction["fluorinated_compounds"] = [val[1] for val in split_values]

        # Average predictions for identical environments
        temp_dataset = {"Code": ["temp"], "SMILES": [smiles]}
        for i in range(71):
            temp_dataset[i] = None
        temp_dataset = pd.DataFrame(temp_dataset, index=[0])

        HOSE_codes = getHoseCodeContent(dataset)
        combined_temp = combined_prediction.merge(HOSE_codes.drop("NMR_Peaks", axis=1), left_index=True, right_index=True)
        combined_temp_grouped = combined_temp.groupby([0, 1, 2, 3, 4, 5])["ensembled_model"].transform("mean")
        combined_temp["ensembled_model"] = combined_temp_grouped
        combined = combined_temp.drop([0, 1, 2, 3, 4, 5], axis=1)
        # Calculate the error of the ensembled model
        combined["ensembled_model_error"] = np.abs(combined["ensembled_model"] - combined["actual"])

        # Usage
        plot_data, table_data, structure_image_base64 = display_results(combined)

        return plot_data, table_data, structure_image_base64

    except Exception as e:
        print(f"Error in predictor: {str(e)}")
        return None, None, None


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_smiles = sys.argv[1]  # Take SMILES string from the command-line argument
        try:
            predictor(input_smiles)  # Attempt to process the SMILES string
        except Exception as e:
            print(f"An error occurred while processing the SMILES string: {e}")
    elif len(sys.argv) == 0:
        smiles = "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
        predictor(smiles)
    else:
        raise ValueError("Please provide a SMILES string as a command-line argument.")
