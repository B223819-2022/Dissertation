"""
This script reads molecules from an SDF file, splits them into individual files, and converts them to PDBQT format using Open Babel. The produced files are placed in a specified directory with each output file named based on the molecule’s ID number.
"""

import os
from openbabel import pybel

def convert_sdf_to_pdbqt(sdf_file, output_dir):
    """
    Convert molecules from an SDF file to PDBQT format and save to an output directory.

    Parameters:
    - sdf_file (str): Path to the input SDF file.
    - output_dir (str): Path to the output directory.

    Side Effects:
    - Creates the output directory if it doesn't exist.
    - Saves converted molecules in the output directory.
    """
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Initialize a counter for the number of molecules processed
    num_molecules = 0

    # Loop through all molecules in the SDF file
    for mol in pybel.readfile('sdf', sdf_file):
        # Get the molecule's ID number
        mol_id = mol.data['IDNUMBER']

        # Use Open Babel to convert the molecule to PDBQT format
        try:
            mol.write('pdbqt', os.path.join(output_dir, f"{mol_id}.pdbqt"), overwrite=True)
            num_molecules += 1
        except Exception as e:
            print(f"Warning: could not process molecule with ID {mol_id}. Error message: {e}")

    # Print the total number of molecules processed
    print(f"Processed {num_molecules} molecules.")

if __name__ == '__main__':
    # Example of how to call the function
    sdf_file = "3D_Library.sdf"
    output_dir = "SD_Library"
    convert_sdf_to_pdbqt(sdf_file, output_dir)
