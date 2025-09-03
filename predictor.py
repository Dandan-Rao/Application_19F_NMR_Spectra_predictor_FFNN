# +
import os
import io
import sys
import base64
import pandas as pd
import pickle

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
    Create .sdf file and temporarily save as /tmp/temp/temp.sdf
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

    # Save .sdf file (use /tmp for Lambda compatibility)
    filename = os.path.join("/tmp", "temp", "temp.sdf")

    w = Chem.SDWriter(filename)
    w.write(sdf)
    w.close()

    # Create enhanced molecular structure image
    file_path = os.path.join("/tmp", "temp", "temp.png")
    
    # Set atom map numbers for display
    for _, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    
    # Create a more appealing molecular structure image with smaller size
    # Use a more reliable highlighting method
    highlight_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "F"]
    highlight_bonds = []
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetSymbol() == "F" or end_atom.GetSymbol() == "F":
            highlight_bonds.append(bond.GetIdx())
    
    # Create color dictionaries with explicit RGB values
    atom_colors = {}
    bond_colors = {}
    
    for atom_idx in highlight_atoms:
        atom_colors[atom_idx] = (1.0, 0.8, 0.0)  # Gold color
    
    for bond_idx in highlight_bonds:
        bond_colors[bond_idx] = (1.0, 0.8, 0.0)  # Gold color
    
    Draw.MolToFile(
        mol, 
        file_path, 
        size=(500, 300),  # Reduced height from 400 to 300 for more compact display
        highlightAtoms=highlight_atoms,
        highlightAtomColors=atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=bond_colors
    )

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
        # Set modern style for better visualization
        plt.style.use("default")
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['axes.edgecolor'] = '#333333'
        
        # Create figure with better proportions and styling
        fig = plt.figure(figsize=(16, 8), facecolor='white')  # Reduced height from 10 to 8
        
        # Create grid layout for better organization (2 rows instead of 3)
        gs = fig.add_gridspec(2, 3, height_ratios=[4, 1], width_ratios=[2, 1, 1], 
                             hspace=0.3, wspace=0.25)
        
        # Main NMR spectrum plot
        ax1 = fig.add_subplot(gs[0, :2])
        
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

        # Enhanced colors and styling
        actual_color = "#2E8B57"  # Sea Green
        predict_color = "#4169E1"  # Royal Blue
        confidence_colors = {
            6: "#00FF00",  # Bright Green
            5: "#32CD32",  # Lime Green
            4: "#FFD700",  # Gold
            3: "#FFA500",  # Orange
            2: "#FF6347",  # Tomato
            1: "#DC143C"   # Crimson
        }

        # Get actual and predicted value counts
        actual = real_PFAS_spectra_df.loc[code, :].value_counts()
        predict = predicted_PFAS_spectra_df.loc[code, :].value_counts()

        # Enhanced NMR spectrum visualization
        if not actual.empty:
            # Plot actual values (negative side)
            ax1.vlines(
                actual.index,
                ymin=0,
                ymax=-actual.values,
                color=actual_color,
                label="Reported Values",
                linewidth=3,
                alpha=0.8,
                capstyle='round'
            )
            # Remove actual value labels - only show peaks

        # Plot predicted values (positive side)
        ax1.vlines(
            predict.index,
            ymin=0,
            ymax=predict.values,
            color=predict_color,
            label="Predicted Values",
            linewidth=3,
            alpha=0.8,
            capstyle='round'
        )
        
        # Add confidence level indicators with colors (only L1-L6 labels)
        for i, j in zip(predict.index, predict.values):
            similarity_level = results[results["ensembled_model"] == i]["similarity_levels"]
            if not similarity_level.empty:
                level = similarity_level.iloc[0]
                color = confidence_colors.get(level, "#808080")
                
                # Enhanced confidence level display (only L1-L6)
                ax1.text(
                    i,
                    j + 0.1,
                    f"L{level}",
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    weight='bold',
                    color='white',
                    bbox=dict(
                        facecolor=color,
                        edgecolor='white',
                        alpha=0.9,
                        boxstyle='round,pad=0.3',
                        linewidth=1
                    )
                )
                
                # Remove predicted value labels - only show L1-L6

        # Set plot limits and styling
        if not actual.empty:
            x_min = min(actual.index.min(), predict.index.min()) - 15
            x_max = max(actual.index.max(), predict.index.max()) + 15
            y_max = max(actual.values.max(), predict.values.max()) + 0.5
        else:
            x_min = predict.index.min() - 15
            x_max = predict.index.max() + 15
            y_max = predict.values.max() + 0.5
            ax1.text(
                (x_min + x_max) / 2,
                -y_max / 2,
                "No reported values available",
                color=actual_color,
                ha="center",
                va="center",
                fontsize=16,
                weight='bold',
                style='italic'
            )

        ax1.set_xlim([x_min, x_max])
        ax1.set_ylim([-y_max, y_max])

        # Enhanced axis labels and styling
        ax1.set_xlabel(r"$^{19}$F NMR Chemical Shift (ppm)", fontsize=14, weight='bold', color='#333333')
        ax1.set_ylabel("Signal Intensity", fontsize=14, weight='bold', color='#333333')
        # Remove title for cleaner appearance
        
        # Remove y-axis ticks and add center line
        ax1.set_yticks([])
        ax1.axhline(0, color='#333333', linewidth=2, alpha=0.7)
        
        # Remove legend for cleaner appearance
        
        # Grid styling
        ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
        ax1.set_facecolor('#F8F9FA')
        
        # Confidence level explanation
        ax2 = fig.add_subplot(gs[0, 2])
        ax2.axis("off")
        
        # Enhanced confidence level table with summary statistics
        confidence_data = {
            "Level": [6, 5, 4, 3, 2, 1],
            "Error": [1.0, 1.8, 4.9, 7.6, 13.7, 13.8],
            "Color": ["#00FF00", "#32CD32", "#FFD700", "#FFA500", "#FF6347", "#DC143C"]
        }
        
        # Calculate summary statistics
        total_predictions = len(results)
        avg_confidence = results['similarity_levels'].mean()
        avg_error = results['ensembled_model_error'].mean()
        
        # Create a clearer and better organized confidence level display with summary
        table_text = "Note: L1-L6 labels above peaks\n"
        table_text += "show confidence levels\n\n"
        table_text += "POSSIBLE ERROR RANGE:\n"
        table_text += "========================\n"
        table_text += "L6: Very High (±1.0 ppm)\n"
        table_text += "L5: High (±1.8 ppm)\n"
        table_text += "L4: Good (±4.9 ppm)\n"
        table_text += "L3: Moderate (±7.6 ppm)\n"
        table_text += "L2: Low (±13.7 ppm)\n"
        table_text += "L1: Very Low (±13.8 ppm)\n\n"
        table_text += "SUMMARY STATISTICS:\n"
        table_text += "========================\n"
        table_text += f"Total Predictions: {total_predictions}\n"
        table_text += f"Avg Confidence: {avg_confidence:.1f}/6\n"
        table_text += f"Avg Error: ±{avg_error:.1f} ppm\n"
        # table_text += f"Compound Code in Dataset: {code if code != 'temp' else 'Not in Dataset'}"
        
        ax2.text(
            0.00, 0.95, table_text, 
            ha="left", va="top",
            fontsize=13,  # Increased from 9 to 11 for better readability
            fontfamily='monospace',
            weight='bold',
            bbox=dict(
                facecolor='#F8F9FA',
                edgecolor='none',  # Removed border by setting edgecolor to 'none'
                alpha=0.95,
                boxstyle='round,pad=0.8',
                linewidth=0  # Set linewidth to 0 to ensure no border
            ),
            transform=ax2.transAxes
        )
        
        # Remove the separate summary statistics section
        # ax3 = fig.add_subplot(gs[1, :]) - This section is now removed

        # Save the enhanced plot
        plot_img = io.BytesIO()
        plt.savefig(plot_img, format="png", dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close(fig)
        plot_img.seek(0)
        plot_base64 = base64.b64encode(plot_img.getvalue()).decode("utf-8")

        # Enhanced table styling
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
                "actual": "Reported Values",
                "ensembled_model": "Predicted Values",
                "ensembled_model_error": "Prediction Error",
            }
        )

        details.reset_index(drop=True, inplace=True)
        
        # Add color coding to confidence levels in the table
        def color_confidence(val):
            if pd.isna(val):
                return ''
            level = int(val)
            colors = {
                6: 'background-color: #90EE90',  # Light Green
                5: 'background-color: #98FB98',  # Pale Green
                4: 'background-color: #FFFFE0',  # Light Yellow
                3: 'background-color: #FFE4B5',  # Moccasin
                2: 'background-color: #FFB6C1',  # Light Pink
                1: 'background-color: #FFC0CB'   # Pink
            }
            return colors.get(level, '')

        # Apply styling to the table
        styled_table = details.style.applymap(color_confidence, subset=['Confidence Level'])
        styled_table = styled_table.set_properties(**{
            'background-color': 'white',
            'color': '#333333',
            'border': '1px solid #ddd',
            'padding': '8px',
            'text-align': 'center'
        }).set_table_styles([
            {'selector': 'th', 'props': [
                ('background-color', '#2C3E50'),
                ('color', 'white'),
                ('font-weight', 'bold'),
                ('text-align', 'center'),
                ('padding', '12px'),
                ('border', '1px solid #34495E')
            ]},
            {'selector': 'td', 'props': [
                ('padding', '10px'),
                ('border', '1px solid #ddd')
            ]}
        ])

        # Save the enhanced table
        table_data = styled_table.to_html(index=False, escape=False)

        # Check if the molecular structure image exists and encode it
        structure_image_path = os.path.join("/tmp", "temp", "temp.png")
        if os.path.exists(structure_image_path):
            with open(structure_image_path, "rb") as structure_file:
                structure_base64 = base64.b64encode(structure_file.read()).decode("utf-8")
        else:
            structure_base64 = None

    except Exception as e:
        print(f"Error generating enhanced results: {str(e)}")
        return None, None, None

    return plot_base64, table_data, structure_base64


# -


def predictor(
    smiles,
    # experimental_data_file_name, 
    train_fluorinated_compounds_file_path=None,
    FFNN_2D_model_path=None, 
    scaler_path=None,
    imputer_path=None, 
    columns_path=None,
    HOSE_Code_database_file_path=None,
):
    # Auto-detect artifact paths for Lambda compatibility
    if train_fluorinated_compounds_file_path is None:
        # Try different possible paths
        possible_paths = [
            os.path.join("artifacts", "Processed_PFAS_19F_NMR_spectra_data.csv"),
            os.path.join("/var/task", "artifacts", "Processed_PFAS_19F_NMR_spectra_data.csv"),
            os.path.join("./artifacts", "Processed_PFAS_19F_NMR_spectra_data.csv")
        ]
        for path in possible_paths:
            if os.path.exists(path):
                train_fluorinated_compounds_file_path = path
                break
        if train_fluorinated_compounds_file_path is None:
            raise FileNotFoundError("Could not find Processed_PFAS_19F_NMR_spectra_data.csv in any expected location")
    
    if FFNN_2D_model_path is None:
        possible_paths = [
            os.path.join("artifacts", "Application_2D_FFNN_sphere5.h5"),
            os.path.join("/var/task", "artifacts", "Application_2D_FFNN_sphere5.h5"),
            os.path.join("./artifacts", "Application_2D_FFNN_sphere5.h5")
        ]
        for path in possible_paths:
            if os.path.exists(path):
                FFNN_2D_model_path = path
                break
        if FFNN_2D_model_path is None:
            raise FileNotFoundError("Could not find Application_2D_FFNN_sphere5.h5 in any expected location")
    
    if scaler_path is None:
        possible_paths = [
            os.path.join("artifacts", "Application_2D_FFNN_scaler_2d_sphere5.pkl"),
            os.path.join("/var/task", "artifacts", "Application_2D_FFNN_scaler_2d_sphere5.pkl"),
            os.path.join("./artifacts", "Application_2D_FFNN_scaler_2d_sphere5.pkl")
        ]
        for path in possible_paths:
            if os.path.exists(path):
                scaler_path = path
                break
        if scaler_path is None:
            raise FileNotFoundError("Could not find Application_2D_FFNN_scaler_2d_sphere5.pkl in any expected location")
    
    if imputer_path is None:
        possible_paths = [
            os.path.join("artifacts", "Application_2D_FFNN_imputer_2d_sphere5.pkl"),
            os.path.join("/var/task", "artifacts", "Application_2D_FFNN_imputer_2d_sphere5.pkl"),
            os.path.join("./artifacts", "Application_2D_FFNN_imputer_2d_sphere5.pkl")
        ]
        for path in possible_paths:
            if os.path.exists(path):
                imputer_path = path
                break
        if imputer_path is None:
            raise FileNotFoundError("Could not find Application_2D_FFNN_imputer_2d_sphere5.pkl in any expected location")
    
    if columns_path is None:
        possible_paths = [
            os.path.join("artifacts", "Application_2D_FFNN_column_names_2d_sphere5.pkl"),
            os.path.join("/var/task", "artifacts", "Application_2D_FFNN_column_names_2d_sphere5.pkl"),
            os.path.join("./artifacts", "Application_2D_FFNN_column_names_2d_sphere5.pkl")
        ]
        for path in possible_paths:
            if os.path.exists(path):
                columns_path = path
                break
        if columns_path is None:
            raise FileNotFoundError("Could not find Application_2D_FFNN_column_names_2d_sphere5.pkl in any expected location")
    
    if HOSE_Code_database_file_path is None:
        possible_paths = [
            os.path.join("artifacts", "HOSE_database_all_fluorianted_compounds.csv"),
            os.path.join("/var/task", "artifacts", "HOSE_database_all_fluorianted_compounds.csv"),
            os.path.join("./artifacts", "HOSE_database_all_fluorianted_compounds.csv")
        ]
        for path in possible_paths:
            if os.path.exists(path):
                HOSE_Code_database_file_path = path
                break
        if HOSE_Code_database_file_path is None:
            raise FileNotFoundError("Could not find HOSE_database_all_fluorianted_compounds.csv in any expected location")
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
