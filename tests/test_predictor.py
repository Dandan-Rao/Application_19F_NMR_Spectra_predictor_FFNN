import pytest
from predictor import *

file_smiles_dict = {
    "43": "C(CC(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)C(=O)O",
    "136": "C(=O)(C(OC(F)(F)F)(F)F)[O-].[Na+]",
    "2": "C(=O)(C(C(F)(F)F)(F)F)O",
    "132": "C(=O)(C(C(F)(F)F)(OC(C(C(F)(F)F)(OC(C(C(F)(F)F)(F)F)(F)F)F)(F)F)F)O",
    "66": "C(C(C(C(F)(F)S(=O)(=O)[O-])(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F.[K+]",
    "59": "C(C(C(F)(F)F)C(F)(F)F)C(=O)O",
    "49": "C(C(=O)O)(C(F)(F)F)F",
    "35": "C(C(=O)O)C(C(F)(F)F)(F)F",
    "63": "C(F)(F)(F)S(=O)(=O)[O-].[Na+]"
}

smiles = []

for _, smile in file_smiles_dict.items():
    smiles.append(smile)

@pytest.mark.parametrize('smile', smiles)
def test_predictor(smile):
    """
    Tests the predictor function with a given SMILES string.

    Args:
        smiles: The SMILES string to test.

    Returns:
        True if the test passes, False otherwise.
    """
    try:
        plot_img, table_data, structure_file = predictor(smile)

        # Assert that the returned data has the expected structure (strings)
        assert isinstance(plot_img, str), f"Expected plot_img to be a string, got {type(plot_img)}"
        assert isinstance(table_data, str), f"Expected table_data to be a string, got {type(table_data)}"
        assert isinstance(structure_file, str), f"Expected structure_file to be a string, got {type(structure_file)}"


        # Perform basic checks on the data
        assert plot_img.strip(), "plot_img is empty after stripping"
        assert table_data.strip(), "table_data is empty after stripping"
        assert structure_file.strip(), "structure_file is empty after stripping"

    except Exception as e:
        pytest.fail(f"Test failed for SMILES '{smile}': {e}")