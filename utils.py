# +
import sys
import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import math
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.impute import KNNImputer
from xgboost import XGBRegressor
from rdkit import Chem
from rdkit.Chem import AllChem

random_state = 1

src_path = os.path.join("external")
if src_path not in sys.path:
    sys.path.append(src_path)


import pickle

sys.path.append(".")
from hosegen import HoseGenerator
from hosegen.geometry import *

# -


# Transform the column names of the DataFrame to integers where possible and keep them as strings otherwise
def convert_column_name(name):
    try:
        return int(name)
    except ValueError:
        return name


# +
def canonical_smiles(smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    return smiles


def Transform_SMILE_to_3D_conformation(SMILES):
    mol = Chem.MolFromSmiles(SMILES)
    hmol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(
        hmol, randomSeed=0xF00D
    )  # attempts to generate a reasonable 3D structure for the molecule by placing the atoms in 3D space, considering the bond lengths, angles, and steric effects.
    AllChem.MMFFOptimizeMolecule(
        hmol
    )  # optimizes the 3D conformation of the molecule using a force field.
    my_embedded_mol = Chem.RemoveHs(hmol)
    return my_embedded_mol


# -

file_path1 = os.path.join("dataset", "descriptors")
file_path2 = os.path.join("dataset", "neighbors")


def Combine_descriptors(
    dataset,
    neighbor_num,
    fluorinated_compounds_CDK_Desc_file_path=file_path1,
    fluorinated_compounds_neighbor_atoms_file_path=file_path2,
    with_additional_info=True,
):
    """
    Parameters
    ----------
    dataset : DataFrame
        Must contain a column 'Code', which will be used to retrieve CDK descriptor files and neighbor atom information from corresponding files.

    neighbor_num : int
        The number of neighboring atoms in 3D space used to extract atomic features.

    fluorinated_compounds_CDK_Desc_file_path : str
        Path to the folder storing CDK descriptor files for all fluorinated compounds. These descriptors are obtained using the CDK library in Java.

    fluorinated_compounds_neighbor_atoms_file_path : str
        Path to the folder storing information about neighboring atoms for all atoms in the fluorinated compounds. This information is obtained using the CDK library in Java.

    with_additional_info: bool
        If True, the 'dataset' DataFrame must also include columns: 'SMILES', 'Compound name', and 'Solvent_used_for_NMR'.

    Returns
    -------
    All_NMR_w_Desc : DataFrame
        A DataFrame containing descriptors from neighboring atoms for each F atom in the fluorinated compounds in the dataset.
    """

    dataset.columns = [convert_column_name(name) for name in dataset.columns]

    fluorinated_compounds_collection = []
    All_NMR_w_Desc = pd.DataFrame()

    # The code will be used to match CDK descriptor files and neighbor atoms file.
    for i in dataset["Code"]:
        fluorinated_compounds_collection.append(i)

    for fluorinated_compounds in fluorinated_compounds_collection:
        # Import Descriptors
        fluorinated_compounds_CDK_Desc_path = os.path.join(
            fluorinated_compounds_CDK_Desc_file_path,
            f"{fluorinated_compounds}_Descriptors.csv",
        )
        fluorinated_compounds_CDK_Desc = pd.read_csv(
            fluorinated_compounds_CDK_Desc_path
        )

        # Delete all columns containg 'proton'
        fluorinated_compounds_CDK_Desc = fluorinated_compounds_CDK_Desc.loc[
            :,
            ~fluorinated_compounds_CDK_Desc.columns.str.contains("proton", case=False),
        ]

        # Get index of all F atoms in the fluorinated_compounds molecular
        F_index = fluorinated_compounds_CDK_Desc[
            fluorinated_compounds_CDK_Desc.SMILES == "F"
        ].index

        Desc_of_target_F_atoms = fluorinated_compounds_CDK_Desc[
            fluorinated_compounds_CDK_Desc.SMILES == "F"
        ]

        # Import information about neighbor atoms
        fluorinated_compounds_neighbor_atom_path = os.path.join(
            fluorinated_compounds_neighbor_atoms_file_path,
            f"{fluorinated_compounds}_Neighbors.csv",
        )
        fluorinated_compounds_neighbor_atoms = pd.read_csv(
            fluorinated_compounds_neighbor_atom_path, header=None
        )

        # Generate a df about NMR peak data, which with F index as index and "NMR_Peaks" as col.
        # If the df does not contains NMR peak data, create a blank df
        try:
            NMR_Peaks = dataset[dataset.Code == fluorinated_compounds]
            NMR_Peaks = NMR_Peaks.reset_index(drop=True)
            NMR_Peaks.index = ["NMR_Peaks"]
            NMR_Peaks = NMR_Peaks.loc[:, F_index]
            NMR_Peaks = NMR_Peaks.T
        except KeyError:
            NMR_Peaks = pd.DataFrame({}, index=F_index)

        # Generate a df about descriptors of target F atoms
        Desc_F_atoms_self = Desc_of_target_F_atoms.rename(columns=lambda x: x + "_self")
        F_NMR_and_Desc = NMR_Peaks.merge(
            Desc_F_atoms_self, left_index=True, right_index=True
        )

        Neighbor_Desc = pd.DataFrame()

        for i in F_index:
            all_neighbors = fluorinated_compounds_neighbor_atoms.shape[1]
            if neighbor_num <= all_neighbors:
                neighbor_index = fluorinated_compounds_neighbor_atoms.iloc[
                    i, range(neighbor_num)
                ]
            else:
                neighbor_index = fluorinated_compounds_neighbor_atoms.iloc[i, :]

            neighbor_desc = fluorinated_compounds_CDK_Desc.iloc[
                neighbor_index
            ].reset_index()
            df = pd.DataFrame()
            for index, _ in neighbor_desc.iterrows():
                interm = pd.DataFrame(neighbor_desc.iloc[index]).T.reset_index()
                interm = interm.rename(columns=lambda x: x + f"_neighbor{index+1}")
                df = df.merge(interm, how="outer", left_index=True, right_index=True)
            df.index = [i]
            Neighbor_Desc = pd.concat([Neighbor_Desc, df], axis=0)
        NMR_w_Desc = pd.concat([F_NMR_and_Desc, Neighbor_Desc], axis=1)
        NMR_w_Desc = NMR_w_Desc.rename(lambda x: f"{x}_{fluorinated_compounds}")

        if with_additional_info:
            # Extract single values for the 'SMILES', 'Compound name', and 'Solvent_used_for_NMR' columns
            NMR_w_Desc["SMILES"] = dataset[dataset.Code == fluorinated_compounds][
                "SMILES"
            ].values[0]
            NMR_w_Desc["Compound name"] = dataset[
                dataset.Code == fluorinated_compounds
            ]["Compound name"].values[0]
            NMR_w_Desc["Solvent_used_for_NMR"] = dataset[
                dataset.Code == fluorinated_compounds
            ]["Solvent_used_for_NMR"].values[0]

        All_NMR_w_Desc = pd.concat([All_NMR_w_Desc, NMR_w_Desc], axis=0)

    return All_NMR_w_Desc


def Fatoms_Peak(
    fluorinated_compounds_NMR_Peaks, fluorinated_compounds_CDK_Desc_file_path=file_path1
):
    fluorinated_compounds_collection = []
    all_NMR_Peaks = pd.DataFrame()
    for i in fluorinated_compounds_NMR_Peaks["Code"]:
        fluorinated_compounds_collection.append(i)

    for fluorinated_compounds in fluorinated_compounds_collection:
        # Import Descriptors
        fluorinated_compounds_CDK_Desc_path = os.path.join(
            fluorinated_compounds_CDK_Desc_file_path,
            f"{fluorinated_compounds}_Descriptors.csv",
        )
        fluorinated_compounds_CDK_Desc = pd.read_csv(
            fluorinated_compounds_CDK_Desc_path
        )

        # Delete all columns containg 'proton'
        fluorinated_compounds_CDK_Desc = fluorinated_compounds_CDK_Desc.loc[
            :,
            ~fluorinated_compounds_CDK_Desc.columns.str.contains("proton", case=False),
        ]

        # Get index of all F atoms in the fluorinated_compounds molecular
        F_index = fluorinated_compounds_CDK_Desc[
            fluorinated_compounds_CDK_Desc.SMILES == "F"
        ].index

        # Generate a df about NMR peak data, which with F index as index and "NMR_Peaks" as col.
        NMR_Peaks = fluorinated_compounds_NMR_Peaks[
            fluorinated_compounds_NMR_Peaks.Code == fluorinated_compounds
        ]
        NMR_Peaks = NMR_Peaks.reset_index(drop=True)
        NMR_Peaks.index = ["NMR_Peaks"]
        NMR_Peaks = NMR_Peaks.loc[:, F_index]
        NMR_Peaks = NMR_Peaks.T

        NMR_Peaks = NMR_Peaks.rename(lambda x: f"{x}_{fluorinated_compounds}")

        # Extract single values for the 'SMILES', 'Compound name', and 'Solvent_used_for_NMR' columns
        NMR_Peaks["SMILES"] = fluorinated_compounds_NMR_Peaks[
            fluorinated_compounds_NMR_Peaks.Code == fluorinated_compounds
        ]["SMILES"].values[0]
        NMR_Peaks["Compound name"] = fluorinated_compounds_NMR_Peaks[
            fluorinated_compounds_NMR_Peaks.Code == fluorinated_compounds
        ]["Compound name"].values[0]
        NMR_Peaks["Solvent_used_for_NMR"] = fluorinated_compounds_NMR_Peaks[
            fluorinated_compounds_NMR_Peaks.Code == fluorinated_compounds
        ]["Solvent_used_for_NMR"].values[0]

        all_NMR_Peaks = pd.concat([all_NMR_Peaks, NMR_Peaks], axis=0)

    return all_NMR_Peaks


# One Hot Encoding of Categorical Columns
def one_hot_encoding_of_smiles_neighborN(
    dataset,
    smiles_neighbors=[
        "SMILES_neighbor5",
        "SMILES_neighbor4",
        "SMILES_neighbor3",
        "SMILES_neighbor2",
        "SMILES_neighbor1",
    ],
):
    encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    cols = smiles_neighbors
    # Fit the encoder to the categorical column(s)
    encoder.fit(dataset[cols])
    encoded = encoder.transform(dataset[cols])

    # Create a DataFrame with the encoded data
    df_encoded = pd.DataFrame(encoded, columns=encoder.get_feature_names_out(cols))
    df_encoded.index = dataset.index

    temp = dataset.drop(columns=cols)

    # Combine the encoded data with the rest of the DataFrame
    dataset_encoded = pd.concat([temp, df_encoded], axis=1)

    print(f"-----Shape of the dataset after encoding the categorical values------")
    print(dataset_encoded.shape)
    return dataset_encoded, encoder


def drop_smiles_neighborN(
    dataset,
    smiles_neighbors=[
        "SMILES_neighbor5",
        "SMILES_neighbor4",
        "SMILES_neighbor3",
        "SMILES_neighbor2",
        "SMILES_neighbor1",
    ],
):
    """
    Drop the categorical values in columns smiles_neighbors
    """
    encoded_columns = [
        col for col in dataset.columns if any(orig in col for orig in smiles_neighbors)
    ]

    return dataset.drop(columns=encoded_columns)


# Data processing
# Drop constant columns
def drop_constant_col(df):
    df = df.loc[:, df.nunique() > 1]
    return df


# Identify bool columns
def find_bool_columns(df):
    bool_columns = df.select_dtypes(exclude=["number"]).sum()
    print("Non-Numeric Columns:")
    print(bool_columns)


# Count non-numeric values in each columns
def count_non_numeric_values(df):
    non_numeric_counts = df.applymap(lambda x: not isinstance(x, (float, int))).sum()

    # Filter to show only columns with non-numeric values more than 0
    filtered_non_numeric_counts = non_numeric_counts[non_numeric_counts > 0]

    print("Columns with Non-Numeric Values:")
    print(filtered_non_numeric_counts)


# Replace all string to NaN
def convert_to_numeric(value):
    try:
        return pd.to_numeric(value)
    except:
        return np.nan


def drop_low_cv_cols(df):
    cv = {}
    for col in df.columns:
        series = df[col]
        std = np.std(series, ddof=1)
        mean = np.mean(series)
        cv[col] = std / mean * 100

    cv_df = pd.DataFrame.from_dict(cv, orient="index")

    print(f"Number of columns with <2% cv: {cv_df[cv_df[0] < 2].shape[0]}")

    cv_df[0].hist(bins=20, range=(0, 100))
    plt.xlabel("CV (%)")
    plt.ylabel("Number of Columns")

    ## Let's drop columns with cv < 2%
    low_cv_cols = cv_df[cv_df[0] < 2].index
    df_new = df.drop(low_cv_cols, axis=1)
    return df_new


def drop_high_ratio_NaN_cols(df):
    bar = 0.8 * df.shape[0]
    df_cleaned = df.loc[:, df.isnull().sum() < bar]
    return df_cleaned


def fill_NaN(X):
    imputer = KNNImputer(n_neighbors=2)
    X_filled = imputer.fit_transform(X)
    X_filled = pd.DataFrame(X_filled, columns=X.columns, index=X.index)
    return X_filled, imputer


def find_highly_correlated_features(df, threshold=0.99):
    correlated_cols = {}
    corr = df.corr()

    for col in df.columns:
        correlated_features = corr.index[corr[col] > threshold].tolist()
        correlated_features.remove(col)  # Remove the feature itself from the list
        if correlated_features:
            correlated_cols[col] = correlated_features
    return correlated_cols


def delete_highly_correlated_features(df, dictionary):
    all_col = df.columns

    for col in all_col:
        if col in df.columns and col in dictionary:
            cols_to_delete = dictionary[col]
            df = df.drop(cols_to_delete, axis=1, errors="ignore")

    return df


def drop_categorical_columns(df):
    # Identify categorical columns (either of type 'object' or 'category')
    categorical_cols = df.select_dtypes(include=["object", "category"]).columns.tolist()

    # Drop the categorical columns from the DataFrame
    df_no_categorical = df.drop(columns=categorical_cols)

    return df_no_categorical


def get_results_table(best_model, X, y):
    y_predict = best_model.predict(X)
    compare = pd.DataFrame(y)
    compare.columns = ["actual"]
    compare["prediction"] = y_predict

    compare["actual"] = pd.to_numeric(compare["actual"])
    compare["prediction"] = pd.to_numeric(compare["prediction"])

    compare["diff"] = compare["actual"] - compare["prediction"]
    compare["diff"] = compare["diff"].abs()
    return compare


def testRidgeCVPerformance(
    dataset, neighbor_num, RidgeCVmodel_path, scaler_path, imputer_path, columns_path
):
    with open(RidgeCVmodel_path, "rb") as file:
        best_model = pickle.load(file)

    with open(scaler_path, "rb") as file:
        scaler = pickle.load(file)

    with open(imputer_path, "rb") as file:
        imputer = pickle.load(file)

    with open(columns_path, "rb") as file:
        train_columns = pickle.load(file)

    # convert values to numeric values when possible
    dataset = Combine_descriptors(
        dataset, neighbor_num=neighbor_num, with_additional_info=True
    )
    #     dataset = dataset.rename_axis('atomCode_fluorinated_compoundsCode', inplace = True)
    dataset.apply(convert_to_numeric)

    # drop rows with NaN values in the 'NMR_Peaks' column
    dataset_dropNaN = dataset.dropna(subset=["NMR_Peaks"])

    dataset_dropNaN = dataset_dropNaN[train_columns]
    dataset_dropNaN_imputed = imputer.transform(dataset_dropNaN)

    dataset_dropNaN_imputed = pd.DataFrame(
        dataset_dropNaN_imputed,
        columns=dataset_dropNaN.columns,
        index=dataset_dropNaN.index,
    )

    y = dataset_dropNaN_imputed["NMR_Peaks"]
    X = dataset_dropNaN_imputed.drop(["NMR_Peaks"], axis=1)

    X_scaled = scaler.transform(X)

    X_scaled = pd.DataFrame(X_scaled)
    X_scaled.columns = X.columns
    X_scaled.index = X.index

    results_table = get_results_table(best_model=best_model, X=X_scaled, y=y)
    plot_prediction_performance(results_table, figure_title=None)
    show_results_scatter(results_table, figure_title=None)
    return results_table


def testRidgePerformance2DFeatures(
    dataset, num_spheres, RidgeCVmodel_path, scaler_path, imputer_path, columns_path
):
    with open(RidgeCVmodel_path, "rb") as file:
        best_model = pickle.load(file)

    with open(scaler_path, "rb") as file:
        scaler = pickle.load(file)

    with open(imputer_path, "rb") as file:
        imputer = pickle.load(file)

    with open(columns_path, "rb") as file:
        train_columns = pickle.load(file)

    get_2d_descriptors = getAtomicDescriptorsFrom2DNeighbors()
    content = get_2d_descriptors.getDescriptorsFromDataset(dataset, num_spheres)

    # Convert all columns in df to numeric where possible, keeping non-numeric values unchanged
    content = content.apply(pd.to_numeric, errors="ignore")

    # Conver column names to 'string'
    content.columns = content.columns.astype(str)

    # Drop rows with NaN values in the 'NMR_Peaks' column
    content = content.dropna(subset=["NMR_Peaks"])

    # Delete columns not shown in the train dataset
    content = content[train_columns]
    content_imputed = imputer.transform(content)

    content_imputed = pd.DataFrame(
        content_imputed, columns=content.columns, index=content.index
    )

    y = content_imputed["NMR_Peaks"]
    X = content_imputed.drop(["NMR_Peaks"], axis=1)

    X_scaled = scaler.transform(X)

    X_scaled = pd.DataFrame(X_scaled)
    X_scaled.columns = X.columns
    X_scaled.index = X.index

    results_table = get_results_table(best_model=best_model, X=X_scaled, y=y)
    plot_prediction_performance(results_table, figure_title=None)
    show_results_scatter(results_table, figure_title=None)
    return results_table


def plot_prediction_performance(results_table, figure_title=None):
    num_below_1, num_below_2, num_below_3, num_below_5, num_below_10, num_below_20 = (
        0,
        0,
        0,
        0,
        0,
        0,
    )
    total_count = len(results_table["diff"])

    for i in results_table["diff"]:
        if i < 20:
            num_below_20 += 1
        if i < 10:
            num_below_10 += 1
        if i < 5:
            num_below_5 += 1
        if i < 3:
            num_below_3 += 1
        if i < 2:
            num_below_2 += 1
        if i < 1:
            num_below_1 += 1

    x = [0, 1, 2, 3, 5, 10, 20]
    y = [
        0,
        num_below_1,
        num_below_2,
        num_below_3,
        num_below_5,
        num_below_10,
        num_below_20,
    ]
    y = [val / total_count for val in y]

    # Create the figure and axes with specified size
    cm = 1 / 2.54  # centimeters in inches
    _, ax = plt.subplots(figsize=(10 * cm, 9 * cm))  # 8 cm x 8 cm

    ax.plot(x, y, marker="o", color="#069AF3")
    ax.set_title(figure_title)
    #     plt.grid(True)
    ax.set_xlim([-1, 22])

    # Set border (edge) line width
    ax.spines["top"].set_linewidth(1)  # Top border
    ax.spines["right"].set_linewidth(1)  # Right border
    ax.spines["bottom"].set_linewidth(1)  # Bottom border
    ax.spines["left"].set_linewidth(1)  # Left border

    # Set axis titles and tick label font sizes
    ax.set_xlabel("Threshold (ppm)", fontsize=14)  # Replace with your label
    ax.set_ylabel("Cumulative Proportion", fontsize=14)  # Replace with your label
    ax.tick_params(axis="x", labelsize=12)  # X-axis numbers font size
    ax.tick_params(axis="y", labelsize=12)  # Y-axis numbers font size

    # Add annotations to each data point
    for i in range(len(x)):
        ax.annotate(
            f"({x[i]}, {y[i]:.2f})",
            (x[i], y[i]),
            textcoords="offset points",
            xytext=(0, 5),
            ha="left",
        )

    plt.show()


def show_results_scatter(results_table, figure_title=None):
    r2 = r2_score(results_table["actual"], results_table["prediction"])

    # Create the figure and axes with specified size
    cm = 1 / 2.54  # centimeters in inches
    _, ax = plt.subplots(figsize=(10 * cm, 9 * cm))

    ax.scatter(
        x=results_table["actual"],
        y=results_table["prediction"],
        alpha=0.6,
        color="#069AF3",
    )

    ax.plot([30, -300], [30, -300], c="red")
    ax.plot([30, -300], [40, -290], c="green", linestyle="dashed")
    ax.plot([30, -300], [20, -310], c="green", linestyle="dashed")

    ax.set_ylim([-280, 30])
    ax.set_xlim([-280, 30])
    ax.set_title(figure_title)
    print(f"R2 = {r2:.2f}")

    # Set border (edge) line width
    ax.spines["top"].set_linewidth(1)  # Top border
    ax.spines["right"].set_linewidth(1)  # Right border
    ax.spines["bottom"].set_linewidth(1)  # Bottom border
    ax.spines["left"].set_linewidth(1)  # Left border

    # Set axis titles and tick label font sizes
    ax.set_xlabel("Actual (ppm)", fontsize=14)  # Replace with your label
    ax.set_ylabel("Prediction (ppm)", fontsize=14)  # Replace with your label
    ax.tick_params(axis="x", labelsize=12)  # X-axis numbers font size
    ax.tick_params(axis="y", labelsize=12)  # Y-axis numbers font size

    mse = mean_squared_error(results_table["actual"], results_table["prediction"])
    rmse = math.sqrt(mse)
    print(f"RMSE = {rmse:.2f}")
    mae = mean_absolute_error(results_table["actual"], results_table["prediction"])
    print(f"MAE = {mae}")
    plt.show()


# +
all_neighbor_atoms_list = [
    "SMILES_neighbor10",
    "SMILES_neighbor9",
    "SMILES_neighbor8",
    "SMILES_neighbor7",
    "SMILES_neighbor6",
    "SMILES_neighbor5",
    "SMILES_neighbor4",
    "SMILES_neighbor3",
    "SMILES_neighbor2",
]


def testModelPerformance_XGBoost_3DNeighborAtoms(
    best_model_file_path, columns_file_path, neighbor_num, test_dataset
):
    """
    Function
    ----------
    Applies the optimized model to the validation or test dataset, outputs the prediction results in a DataFrame,
    and visualizes the prediction performance using both a line plot and a scatter plot.

    Parameters
    ----------
    best_model_file_path: Path to the .json file containing the optimized model.

    columns_file_path: Path to the .pkl file containing the column names of the dataset used for modeling.

    test_dataset: A DataFrame containing information about fluorinated compounds and their corresponding 19F NMR shift values.

    Output
    ----------
    results_table: A DataFrame containing the prediction results for the test_dataset.
    """
    best_model = XGBRegressor()
    best_model.load_model(best_model_file_path)

    with open(columns_file_path, "rb") as f:
        train_cols = pickle.load(f)

    # Step 1. Combine NMR shift values of F atoms with its features
    fluorinated_compounds_w_Desc = Combine_descriptors(
        test_dataset, neighbor_num=neighbor_num
    )

    # Step 2. Only keep columns that were used in the dataset for modeling while delete other columns
    fluorinated_compounds_w_Desc = fluorinated_compounds_w_Desc[train_cols]

    # For validation purpose, drop rows with NaN values in the 'NMR_Peaks' column
    fluorinated_compounds_w_Desc = fluorinated_compounds_w_Desc.dropna(
        subset=["NMR_Peaks"]
    )

    # Get y values
    y = fluorinated_compounds_w_Desc["NMR_Peaks"]

    # Drop non-numeric columns from X
    orig_features = ["NMR_Peaks"]
    X = fluorinated_compounds_w_Desc.drop(orig_features, axis=1)

    # Ensure all values in the X are numerical values
    X = X.apply(pd.to_numeric)

    results_table = get_results_table(best_model=best_model, X=X, y=y)

    plot_prediction_performance(results_table, figure_title=None)

    show_results_scatter(results_table, figure_title=None)

    return results_table


# -


class getAtomicDescriptorsFrom2DNeighbors:
    # Define a function to get a list of atom indexs in the next spheres
    def findNeighborsNextSphere(
        self, last_sphere_atoms, current_sphere_atoms, next_sphere_level
    ):
        """Get a list of atom indexs in the next spheres

        Parameters
        ----------
        last_sphere_atoms: list
            A list of atoms in the last sphere.

        current_sphere_atoms: list
            A list of atoms in the current sphere

        next_sphere_level: int
             The number of next sphere. For target F atom, the C atom next to it belongs to sphere 0.

        Output
        ----------
        A list of indexs of atoms in the next sphere. If the sphere is No.4. Then the output list will be a list of length 12.
        """
        neighbor_list = []

        # Convert last_sphere_atoms to a set of indices for faster lookups
        if last_sphere_atoms is not None:
            last_sphere_atoms_set = set(
                atom.GetIdx() for atom in last_sphere_atoms if atom is not None
            )
        else:
            last_sphere_atoms_set = set()

        for atom in current_sphere_atoms:
            if atom is None:
                neighbors = [
                    None
                ] * 3  # Adjusted to 3 neighbors. In our case, each atom have at most 4 neighbor atoms.
            else:
                neighbors = (
                    atom.GetNeighbors()
                )  # This function will get all neighbor atoms
                # Filter out atoms that are in last_sphere_atoms
                neighbors = [
                    a for a in neighbors if a.GetIdx() not in last_sphere_atoms_set
                ]

                # Ensure the number of neighbor atoms in the next sphere for each atom is exactly 3.
                if len(neighbors) < 3:
                    neighbors.extend([None] * (3 - len(neighbors)))
                elif len(neighbors) > 3:
                    print(f"More than 3 neighbors: {len(neighbors)}. Trimming to 3.")
                    neighbors = neighbors[:3]

            neighbor_list.extend(neighbors)

        # Final length should be 3^sphere_level
        final_length = 3**next_sphere_level
        if len(neighbor_list) < final_length:
            neighbor_list.extend([None] * (final_length - len(neighbor_list)))

        return neighbor_list

    # Get neighbors of each F atom
    def getNeighborsOfFAtom(self, F_atom):
        """For the target F atom, get neighbor atoms in the 6 neighbor spheres
        Parameters
        ----------
        F_atom: atom object

        Output
        ----------
        sphere0[:1], sphere1, sphere2, sphere3, sphere4, sphere5: list, list, list, list, list, list
            index of atoms in sphere 0, index of atoms in sphere 1, index of atoms in sphere 2, index of atoms in sphere 3,
            index of atoms in sphere 4, and index of atoms in sphere 5
        """
        sphere0_atoms = self.findNeighborsNextSphere(None, [F_atom], 0)
        sphere1_atoms = self.findNeighborsNextSphere(
            [F_atom], sphere0_atoms[:1], 1
        )  # We know sphere0 can only have one valid atom, C
        sphere2_atoms = self.findNeighborsNextSphere(
            sphere0_atoms[:1], sphere1_atoms, 2
        )
        sphere3_atoms = self.findNeighborsNextSphere(sphere1_atoms, sphere2_atoms, 3)
        sphere4_atoms = self.findNeighborsNextSphere(sphere2_atoms, sphere3_atoms, 4)
        sphere5_atoms = self.findNeighborsNextSphere(sphere3_atoms, sphere4_atoms, 5)

        sphere0 = [
            atom.GetIdx() if atom is not None else None for atom in sphere0_atoms
        ]
        sphere1 = [
            atom.GetIdx() if atom is not None else None for atom in sphere1_atoms
        ]
        sphere2 = [
            atom.GetIdx() if atom is not None else None for atom in sphere2_atoms
        ]
        sphere3 = [
            atom.GetIdx() if atom is not None else None for atom in sphere3_atoms
        ]
        sphere4 = [
            atom.GetIdx() if atom is not None else None for atom in sphere4_atoms
        ]
        sphere5 = [
            atom.GetIdx() if atom is not None else None for atom in sphere5_atoms
        ]

        return (
            sphere0[:1],
            sphere1,
            sphere2,
            sphere3,
            sphere4,
            sphere5,
        )  # Only keep the only one valid atom, C

    def getNeighborsInDiffSpheres(self, smiles):
        """
        Parameter
        ----------
        smiles: string
            SMILES of the target compound

        Output
        ----------
        df： Dataframe
            Each line in the df shows the index of neighbor atoms for one F atom in the molecule.
            column name of the df are the number of spheres: 0, 1, 2, 3, 4, 5.
            index of the df are index of F atoms
        """
        # Create an RDKit molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        neighbors = {}
        for atom in mol.GetAtoms():
            atom_symbol = (
                atom.GetSymbol()
            )  # Get the atomic symbol (e.g., 'C', 'O', 'F')
            if atom_symbol == "F":
                atom_index = atom.GetIdx()  # Get the index of the atom
                sphere0, sphere1, sphere2, sphere3, sphere4, sphere5 = (
                    self.getNeighborsOfFAtom(atom)
                )
                neighbors[atom_index] = [
                    sphere0,
                    sphere1,
                    sphere2,
                    sphere3,
                    sphere4,
                    sphere5,
                ]

        df = pd.DataFrame.from_dict(neighbors).T
        df = df.applymap(lambda x: tuple(x) if isinstance(x, list) else x)
        df = df.drop_duplicates()
        return df

    def getFeaturesTableForAtoms(
        self,
        smiles,
        feature_list=[
            "mass",
            "hybridization",
            "isAromatic",
            "degree",
            "valence",
            "explicit_valence",
            "isInRing",
        ],
    ):
        """Select some atomic features, then generate a table showing these features of each atom in the target molecule.
        Parameters
        ----------
        smiles: string
            SMILES of the target molecule
        feature_list: list[string]
            Controls which features will be remained

        Output
        ----------
        df: DataFrame
            columns are atomic features
            index are atom index in the molecule
        """
        mol = Chem.MolFromSmiles(smiles)
        atom_prop = {}
        for atom in mol.GetAtoms():
            atom_symbol = atom.GetSymbol()
            atom_index = atom.GetIdx()
            #             atom_atomicNum = atom.GetAtomicNum()
            atom_mass = atom.GetMass()
            atom_hybridization = atom.GetHybridization()  # SP1, SP2, SP3
            atom_isAromatic = atom.GetIsAromatic()  # bool
            atom_degree = (
                atom.GetDegree()
            )  # The degree of an atom is defined to be its number of directly-bonded neighbors.
            atom_valence = (
                atom.GetTotalValence()
            )  # → int. Returns the total valence (explicit + implicit) of the atom.
            atom_explicit_valence = atom.GetExplicitValence()
            atom_isInRing = atom.IsInRing()
            atom_prop[atom_index] = [
                atom_symbol,
                atom_mass,
                atom_hybridization,
                atom_isAromatic,
                atom_degree,
                atom_valence,
                atom_explicit_valence,
                atom_isInRing,
            ]
            # transform to a dataframe
            df = pd.DataFrame.from_dict(atom_prop).T
            df.columns = [
                "symbol",
                "mass",
                "hybridization",
                "isAromatic",
                "degree",
                "valence",
                "explicit_valence",
                "isInRing",
            ]
            # Convert the column to strings
            df["hybridization"] = df["hybridization"].astype(str).str.slice(-1)
            # Convert the last character to integers, using errors='coerce' to handle invalid conversion
            df["hybridization"] = pd.to_numeric(df["hybridization"], errors="coerce")
            df[["isAromatic", "isInRing"]] = df[["isAromatic", "isInRing"]].astype(int)
            df = df[feature_list]
        return df

    def removeNoneFromTuple(self, tup):
        """Filter out None values from the tuple and return a new tuple"""
        return [x for x in tup if x is not None]

    def slimNeighbors(self, neighbors):
        """Apply the function to each element of the DataFrame"""
        return neighbors.applymap(self.removeNoneFromTuple)

    def getNeighborsList(self, slim_neighbors, num_spheres):
        # sphere 0 == 1 atom
        # sphere 1 <= 3 atoms
        # sphere 2, leave 9 spaces
        # sphere 3, leave 10 spaces
        # sphere 4, leave 10 spaces
        # if more atoms appear in sphere 2~5, keep only the first 9 or 10 atoms.
        """Reshape the neighbor atoms dataframe for modeling purpose.
        Parameters
        ----------
        slim_neighbors: DataFrame
            A dataframe with each row shows the list of atom indexs in each neighbor spheres of a specific F atom

        Output
        ----------
        df: DataFrame
            A DataFrame with only 33 columns, each column is a neighbor atom slot
            Index of the DataFrame if the atom indexs of F atoms in the molecule.
        """
        dic = {}
        for _, row in slim_neighbors.iterrows():
            neighbors_list = [None] * (1 + 3 + 9 + 10 + 10)
            neighbors_list[0] = row[0][0]
            # Handle sphere 1: up to 3 atoms
            for j in range(min(len(row[1]), 3)):
                neighbors_list[1 + j] = row[1][j]

            # Handle sphere 2: up to 9 atoms
            for j in range(min(len(row[2]), 9)):
                neighbors_list[4 + j] = row[2][j]

            # Handle sphere 3: up to 10 atoms
            for j in range(min(len(row[3]), 10)):
                neighbors_list[13 + j] = row[3][j]

            # Handle sphere 4: up to 10 atoms
            for j in range(min(len(row[4]), 10)):
                neighbors_list[23 + j] = row[4][j]

            # Store the neighbors list in the dictionary
            dic[i] = neighbors_list
            df = pd.DataFrame.from_dict(dic).T
            df.columns = [
                "0_1",
                "1_1",
                "1_2",
                "1_3",
                "2_1",
                "2_2",
                "2_3",
                "2_4",
                "2_5",
                "2_6",
                "2_7",
                "2_8",
                "2_9",
                "3_1",
                "3_2",
                "3_3",
                "3_4",
                "3_5",
                "3_6",
                "3_7",
                "3_8",
                "3_9",
                "3_10",
                "4_1",
                "4_2",
                "4_3",
                "4_4",
                "4_5",
                "4_6",
                "4_7",
                "4_8",
                "4_9",
                "4_10",
            ]
            if num_spheres == 1:
                df = df[["0_1"]]
            if num_spheres == 2:
                df = df[["0_1", "1_1", "1_2", "1_3"]]
            if num_spheres == 3:
                df = df[
                    [
                        "0_1",
                        "1_1",
                        "1_2",
                        "1_3",
                        "2_1",
                        "2_2",
                        "2_3",
                        "2_4",
                        "2_5",
                        "2_6",
                        "2_7",
                        "2_8",
                        "2_9",
                    ]
                ]
            if num_spheres == 4:
                df = df[
                    [
                        "0_1",
                        "1_1",
                        "1_2",
                        "1_3",
                        "2_1",
                        "2_2",
                        "2_3",
                        "2_4",
                        "2_5",
                        "2_6",
                        "2_7",
                        "2_8",
                        "2_9",
                        "3_1",
                        "3_2",
                        "3_3",
                        "3_4",
                        "3_5",
                        "3_6",
                        "3_7",
                        "3_8",
                        "3_9",
                        "3_10",
                    ]
                ]
            if num_spheres == 5:
                df = df[
                    [
                        "0_1",
                        "1_1",
                        "1_2",
                        "1_3",
                        "2_1",
                        "2_2",
                        "2_3",
                        "2_4",
                        "2_5",
                        "2_6",
                        "2_7",
                        "2_8",
                        "2_9",
                        "3_1",
                        "3_2",
                        "3_3",
                        "3_4",
                        "3_5",
                        "3_6",
                        "3_7",
                        "3_8",
                        "3_9",
                        "3_10",
                        "4_1",
                        "4_2",
                        "4_3",
                        "4_4",
                        "4_5",
                        "4_6",
                        "4_7",
                        "4_8",
                        "4_9",
                        "4_10",
                    ]
                ]
        return df

    def getDescriptorsList(self, neighbors_list, atom_features):
        """Combine the list of neighbor atoms with atomic features, create a new df
        Parameters
        ----------
        smiles: string
            SMILES of the target molecule

        Output
        ----------
        df: DataFrame
            A new dataframe contains neighbor atoms and thier atomic features.
        """
        dic = {}
        for i, row in neighbors_list.iterrows():
            descriptors = []
            for neighbor in row:
                if isinstance(neighbor, (int, float)) and not pd.isna(neighbor):
                    neighbor = int(neighbor)
                    descriptors.extend(atom_features.iloc[neighbor, :].values)
                else:
                    # Append a list of Nones if there is no neighbor (at most 6 features as per the example)
                    descriptors.extend([None] * atom_features.shape[1])

            dic[i] = descriptors
        df = pd.DataFrame.from_dict(dic).T
        return df

    def getDescriptorsFromSmiles(
        self,
        smiles,
        num_spheres,
        feature_list=[
            "mass",
            "hybridization",
            "isAromatic",
            "degree",
            "valence",
            "explicit_valence",
            "isInRing",
        ],
    ):
        """Combine all above functions, finish all steps in one function
        Parameters
        ----------
        smiles: string
            SMILES of the target molecule

        Output
        ----------
        df: DataFrame
            A dataframe contains neighbor atoms and thier atomic features.
        """
        feature_table = self.getFeaturesTableForAtoms(smiles, feature_list)
        neighbors_in_diff_spheres = self.getNeighborsInDiffSpheres(smiles)
        slim_neighbors = self.slimNeighbors(neighbors_in_diff_spheres)
        neighbors_list = self.getNeighborsList(slim_neighbors, num_spheres)
        content = self.getDescriptorsList(neighbors_list, feature_table)
        return content

    def getDescriptorsFromDataset(
        self,
        dataset,
        num_spheres,
        feature_list=[
            "mass",
            "hybridization",
            "isAromatic",
            "degree",
            "valence",
            "explicit_valence",
            "isInRing",
        ],
    ):
        """Get a dataframe with each row being the atomic features of neighbor atoms of a target F atom in a fluorinated compound"""
        dataset.columns = [convert_column_name(name) for name in dataset.columns]
        fluorinated_compounds_content = pd.DataFrame()
        for _, row in dataset.iterrows():
            smiles = row["SMILES"]
            fluorinated_compounds = row["Code"]
            content = self.getDescriptorsFromSmiles(smiles, num_spheres, feature_list)
            index_list = content.index
            try:
                content["NMR_Peaks"] = row[index_list]
            except (KeyError, IndexError):
                pass

            content = content.rename(lambda x: f"{x}_{fluorinated_compounds}")
            fluorinated_compounds_content = pd.concat(
                [fluorinated_compounds_content, content], axis=0
            )
        return fluorinated_compounds_content

    def testXGBoost2DModelPerformance(
        self,
        best_model_file_path,
        dataset,
        num_spheres,
        feature_list=[
            "mass",
            "hybridization",
            "isAromatic",
            "degree",
            "valence",
            "explicit_valence",
            "isInRing",
        ],
    ):
        best_model = XGBRegressor()
        best_model.load_model(best_model_file_path)

        get_2d_descriptors = getAtomicDescriptorsFrom2DNeighbors()
        vali_content = get_2d_descriptors.getDescriptorsFromDataset(
            dataset, num_spheres, feature_list
        )

        vali_content = vali_content.dropna(subset=["NMR_Peaks"])
        y = vali_content["NMR_Peaks"]
        X = vali_content.drop(["NMR_Peaks"], axis=1)

        results_table = get_results_table(best_model=best_model, X=X, y=y)
        plot_prediction_performance(results_table, figure_title=None)
        show_results_scatter(results_table, figure_title=None)
        return results_table


## Generate HOSE code with different radius
def getHoseCodeContent(data, max_radius=6):
    """For each fluorinated compound in the dataset, retrieve the HOSE code for the F atoms in the molecule.
    Parameters
    ----------
    data: DataFrame
        A dataframe contains SMILES, Code of various fluorinated compounds

    Output
    ----------
    fluorinated_compounds_content: DataFrame
        The index of the DataFrame corresponds to the F atom code, formatted as `'FatomIndex_fluorinated_compoundsCode'`.
        The columns represent radii from 0 to 5.
        For example, the values in column 3 contain the HOSE code that encodes information about neighboring atoms within 3 spheres.
        The final column contains the 19F NMR chemical shift values.
    """
    gen = HoseGenerator()
    # Transform the column names of the DataFrame to integers where possible and keep them as strings otherwise
    data.columns = [convert_column_name(name) for name in data.columns]
    fluorinated_compounds_content = pd.DataFrame()

    for _, row in data.iterrows():
        smiles = row["SMILES"]
        fluorinated_compounds = row["Code"]
        mol = Chem.MolFromSmiles(smiles)
        output_file_path = f"..//temp//{fluorinated_compounds}.mol"
        Chem.MolToMolFile(mol, output_file_path)
        wedgemap = create_wedgemap(output_file_path)
        if os.path.exists(output_file_path):
            os.remove(output_file_path)

        dic = {}
        for atom in mol.GetAtoms():
            atom_symbol = atom.GetSymbol()
            atom_index = atom.GetIdx()
            if atom_symbol == "F":
                hose_codes = []
                for j in range(max_radius):
                    hose_codes.append(
                        gen.get_Hose_codes(
                            mol,
                            atom_index,
                            usestereo=True,
                            max_radius=j + 1,
                            wedgebond=wedgemap,
                        )
                    )
                dic[atom_index] = hose_codes

        df = pd.DataFrame.from_dict(dic).T
        F_indexs = df.index
        df["NMR_Peaks"] = row[F_indexs].values
        df = df.rename(index=lambda x: f"{x}_{fluorinated_compounds}")
        #         df = df.drop_duplicates()
        fluorinated_compounds_content = pd.concat(
            [fluorinated_compounds_content, df], axis=0
        )
    return fluorinated_compounds_content


# +
def getTrainDictionary_HOSE(train_data):
    """
    Get the mean value of 'NMR_Peaks' for each unique values in HOSE code.
    Parameters
    ----------
    train_data: DataFrame
        A datafram with F atom index being df index and radius being column-name. Values are HOSE code of different radius.

    Output
    ----------
    [sphere6_dic, sphere5_dic, sphere4_dic, sphere3_dic, sphere2_dic, sphere1_dic]: A list of dictionaries
    """
    train_data = train_data.dropna(subset=["NMR_Peaks"])
    train_data["NMR_Peaks"] = train_data["NMR_Peaks"].astype(float)

    # Get the mean value of 'NMR_Peaks' for each unique values in column '5'
    grouped_df = train_data[[5, "NMR_Peaks"]].groupby(5)["NMR_Peaks"].mean()
    sphere6_dic = grouped_df.to_dict()

    grouped_df = train_data[[4, "NMR_Peaks"]].groupby(4)["NMR_Peaks"].mean()
    sphere5_dic = grouped_df.to_dict()

    grouped_df = train_data[[3, "NMR_Peaks"]].groupby(3)["NMR_Peaks"].mean()
    sphere4_dic = grouped_df.to_dict()

    grouped_df = train_data[[2, "NMR_Peaks"]].groupby(2)["NMR_Peaks"].mean()
    sphere3_dic = grouped_df.to_dict()

    grouped_df = train_data[[1, "NMR_Peaks"]].groupby(1)["NMR_Peaks"].mean()
    sphere2_dic = grouped_df.to_dict()

    grouped_df = train_data[[0, "NMR_Peaks"]].groupby(0)["NMR_Peaks"].mean()
    sphere1_dic = grouped_df.to_dict()

    return [
        sphere6_dic,
        sphere5_dic,
        sphere4_dic,
        sphere3_dic,
        sphere2_dic,
        sphere1_dic,
    ]


def HOSE_Model(sphere_dics, test_data, mean_value_in_train_data):
    """
    Parameters
    ----------
    sphere_dics: A list of dictionaries
        With max_radius = 6, the list contains 6 dictionary.
        The key of the nth dictionary is the HOSE code with radius being n, and values being the mean 19F NMR
        shift value of the HOSE code in the training dataset.

    test_data: DataFrame
        A dataframe contains the HOSE code for each F atoms in the test dataset.

    mean_value_in_train_data: float
        If no same HOSE code was found in the sphere_dics, use mean_value_in_train_data as prediction value

    Output
    ----------
    prediction: list
        The predicted 19F NMR shift values for F atoms in the test dataset

    similarity_levels: list[int]
        Meaning of the similarity level:
            6: find in pool with max_radius = 6
            5: find in pool with max_radius = 5

            # HOSE code prediction with four or more spheres and respecting stereochemistry are generally considered reliable.
            4: find in pool with max_radius = 4
            3: find in pool with max_radius = 3
            2: find in pool with max_radius = 2
            1: find in pool with max_radius = 1
    """
    test_data["NMR_Peaks"] = test_data["NMR_Peaks"].apply(
        pd.to_numeric, errors="coerce"
    )

    prediction = []
    similarity_levels = []
    for _, row in test_data.iterrows():
        if row[5] in sphere_dics[0]:
            prediction.append(sphere_dics[0][row[5]])
            similarity_levels.append(6)
        elif row[4] in sphere_dics[1]:
            prediction.append(sphere_dics[1][row[4]])
            similarity_levels.append(5)
        elif row[3] in sphere_dics[2]:
            prediction.append(sphere_dics[2][row[3]])
            similarity_levels.append(4)
        elif row[2] in sphere_dics[3]:
            prediction.append(sphere_dics[3][row[2]])
            similarity_levels.append(3)
        elif row[1] in sphere_dics[4]:
            prediction.append(sphere_dics[4][row[1]])
            similarity_levels.append(2)
        elif row[0] in sphere_dics[5]:
            prediction.append(sphere_dics[5][row[0]])
            similarity_levels.append(1)
        else:
            prediction.append(mean_value_in_train_data)
            similarity_levels.append(0)
    return prediction, similarity_levels


def getResults_HOSE(prediction, similarity_levels, test_data):
    """Reorganize prediction, similarity_levels and combine it with test_data
    Output
    ----------
    results: DataFrame
        Column 'actual': Actual 19F NMR shift values for F atoms in the test dataset
        Column 'prediction': The predicted 19F NMR shift values using the HOSE code method
        Column 'diff': Absolute  error of the prediction
        Column 'similarity_levels': The value Represents how similar a particular F atom is to the F atoms in the training dataset.
        This can be treated as a value between 1 and 6, where 6 represents maximum similarity (i.e., very reliable predictions).
    """
    test_data["NMR_Peaks"] = test_data["NMR_Peaks"].apply(
        pd.to_numeric, errors="coerce"
    )
    results = test_data[["NMR_Peaks"]].copy()
    results.columns = ["actual"]
    results["prediction"] = prediction
    results["diff"] = results["prediction"] - results["actual"]
    results["diff"] = results["diff"].abs()
    results["similarity_levels"] = similarity_levels
    return results
