"""
This script uses PSOVina to dock a set of ligands (from an input directory) to a specified protein using the coordinates of a specified pocket and places the results in a directory. The script utilises parallels and a default exhaustiveness of 8 and creates individual log files for each molecule docking and places them in a directory.   
The docking results files are then analysed and sorted based on binding affinity, and the top 500 ligands are copied to a designated directory.
"""

import os
import shutil
from subprocess import Popen
import subprocess
from tqdm import tqdm
import psutil

def docking_process(sd_library_dir, protein_file, center_x, center_y, center_z, size_x, size_y, size_z):
    """
    Perform docking of ligands to a protein using PSOVina and sort results by binding affinity.
    
    Parameters:
    - sd_library_dir (str): Directory containing ligand files in PDBQT format.
    - protein_file (str): Path to the protein file in PDBQT format.
    - center_x, center_y, center_z (float): Coordinates for the center of the docking pocket.
    - size_x, size_y, size_z (float): Dimensions of the docking pocket.

    Side Effects:
    - Creates result directories if they don't exist.
    - Saves docking results in the designated directory.
    - Copies the top 500 ligands with the best binding affinities to a separate directory.
    """
    
    # Define the directories and filenames
    results_dir = "dock1_results_2bv6"
    top500_dir = "dock1_top500_2bv6"
    log_dir = "log_dock1_2bv6"

    # Create the directories if they don't exist
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(top500_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Prepare the PSOVina command template
    psovina_cmd_template = "./psovina --receptor {} --ligand {} --center_x {} --center_y {} --center_z {} --size_x {} --size_y {} --size_z {} --log {} --out {} --exhaustiveness 8"

    # Get the list of PDBQT files in the SD_Library directory
    pdbqt_files = [f for f in os.listdir(sd_library_dir) if f.endswith('.pdbqt')]

    # Initialize the progress bar
    pbar = tqdm(total=len(pdbqt_files), desc="Docking progress", mininterval=0.5, ascii=True)

    # Run the PSOVina commands in parallel
    processes = set()
    for pdbqt_file in pdbqt_files:
        results_file = os.path.join(results_dir, pdbqt_file)
        log_file_path = os.path.join(log_dir, f"log_{pdbqt_file[:-6]}.txt")
        psovina_cmd = psovina_cmd_template.format(protein_file, os.path.join(sd_library_dir, pdbqt_file), center_x, center_y, center_z, size_x, size_y, size_z, log_file_path, results_file)
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

    # Collect the results and sort by binding affinity
    results = []
    for pdbqt_file in os.listdir(results_dir):
        if pdbqt_file.endswith('.pdbqt'):
            with open(os.path.join(results_dir, pdbqt_file), 'r') as file:
                for line in file:
                    if line.startswith("REMARK VINA RESULT:"):
                        binding_affinity = float(line.split()[3])
                        results.append((binding_affinity, pdbqt_file))
                        break
    results.sort(key=lambda x: x[0])

    # Copy the original PDBQT files with the top 500 best binding affinities to the top500 directory
    for _, pdbqt_file in results[:500]:
        shutil.copy(os.path.join(sd_library_dir, pdbqt_file), top500_dir)

    print("Docking process completed. Top 500 molecules have been copied to the specified directory.")

if __name__ == '__main__':
    # Example of how to call the function with specific parameters
    sd_library_dir = "SD_Library"
    protein_file = "2bv6_final_proto.pdbqt"
    center_x = 80.64
    center_y = 2.39
    center_z = 4.29
    size_x = 30
    size_y = 30
    size_z = 30
    docking_process(sd_library_dir, protein_file, center_x, center_y, center_z, size_x, size_y, size_z)
