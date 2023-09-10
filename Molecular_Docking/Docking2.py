"""
This script uses PSOVina to perform a second round of docking using the specified directory containing the top 500 ligands from the initial docking. The script utilises parallels and an exhaustiveness of 32 and creates individual log files for the docking of each molecule and places them in a directory.   
"""

import os
import shutil
from subprocess import Popen
import subprocess
from tqdm import tqdm
import psutil

def second_round_docking(input_dir, protein_file, center_x, center_y, center_z, size_x, size_y, size_z):
    """
    Perform a second round of docking of ligands to a protein using PSOVina.
    
    Parameters:
    - input_dir (str): Directory containing ligand files in PDBQT format for the second round.
    - protein_file (str): Path to the protein file in PDBQT format.
    - center_x, center_y, center_z (float): Coordinates for the center of the docking pocket.
    - size_x, size_y, size_z (float): Dimensions of the docking pocket.

    Side Effects:
    - Creates result directories if they don't exist.
    - Saves docking results and log files in the designated directories.
    """
    
    # Define the directories and filenames
    results_dir = "dock2_results_2bv6"
    log_dir = "log_dock2_2bv6"

    # Create the directories if they don't exist
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Prepare the PSOVina command template
    psovina_cmd_template = "./psovina --receptor {} --ligand {} --center_x {} --center_y {} --center_z {} --size_x {} --size_y {} --size_z {} --log {} --out {} --exhaustiveness 32"

    # Get the list of PDBQT files in the input directory
    pdbqt_files = [f for f in os.listdir(input_dir) if f.endswith('.pdbqt')]

    # Initialize the progress bar
    pbar = tqdm(total=len(pdbqt_files), desc="Docking progress", mininterval=0.5, ascii=True)

    # Run the PSOVina commands in parallel
    processes = set()
    for pdbqt_file in pdbqt_files:
        results_file = os.path.join(results_dir, pdbqt_file)
        log_file_path = os.path.join(log_dir, f"log_{pdbqt_file[:-6]}.txt")
        psovina_cmd = psovina_cmd_template.format(protein_file, os.path.join(input_dir, pdbqt_file), center_x, center_y, center_z, size_x, size_y, size_z, log_file_path, results_file)
        processes.add(Popen(psovina_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT))
        # Limit the number of parallel processes to the number of available CPU cores
        if len(processes) >= psutil.cpu_count():
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])
        pbar.update()

    pbar.close()

    # Wait for all processes to finish
    for p in processes:
        if p.poll() is None:
            p.wait()

    print("Second round of docking completed. Results and log files are saved in the specified directories.")

if __name__ == '__main__':
    # Example of how to call the function with specific parameters
    input_dir = "dock1_top500_2bv6"
    protein_file = "2bv6_final_proto.pdbqt"
    center_x = 80.64
    center_y = 2.39
    center_z = 4.29
    size_x = 30
    size_y = 30
    size_z = 30
    second_round_docking(input_dir, protein_file, center_x, center_y, center_z, size_x, size_y, size_z)
