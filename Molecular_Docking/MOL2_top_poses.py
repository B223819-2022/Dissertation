"""
This script analyses the results from two rounds of docking and extracts the top 10 poses and converts them to the MOL2 format using Open Babel (molecules that produce warning messages during the conversion are ignored). The script also compiles a score file that documents the binding affinity scores for each molecule from the first and second dockings. 
"""

import os
import subprocess

def get_best_score_and_pose(file_path):
    """
    Extracts the best binding affinity score and the corresponding pose from a docking result file.
    
    Parameters:
    - file_path (str): Path to the docking result file (in PDBQT format).
    
    Returns:
    - float: Best binding affinity score.
    - list: Lines representing the best pose in the file.
    """
    best_score = None
    best_pose_lines = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if "REMARK VINA RESULT:" in line:
                best_score = float(line.split()[3])
                break
        for line in lines:
            if "MODEL" in line:
                best_pose_lines = lines[lines.index(line):lines.index("ENDMDL\n")+1]
                break
    return best_score, best_pose_lines

def analyze_results(results_dir1, results_dir2, final_poses_dir, scores_file):
    """
    Analyzes docking results, extracts top poses, converts them to MOL2 format, and compiles a scores file.
    
    Parameters:
    - results_dir1 (str): Directory containing results from the first round of docking.
    - results_dir2 (str): Directory containing results from the second round of docking.
    - final_poses_dir (str): Directory where the top poses from the second round will be saved in MOL2 format.
    - scores_file (str): Path to the output file where scores will be documented.
    """
    
    os.makedirs(final_poses_dir, exist_ok=True)

    results = []
    for pdbqt_file in os.listdir(results_dir2):
        if pdbqt_file.endswith('.pdbqt'):
            second_docking_score, best_pose_lines = get_best_score_and_pose(os.path.join(results_dir2, pdbqt_file))
            if second_docking_score is None:
                print(f"Failed to extract binding affinity score from {pdbqt_file}")
            results.append((second_docking_score, pdbqt_file, best_pose_lines))

    results.sort(key=lambda x: x[0])  # Sorting by binding affinity, most negative (best) first
    top_10_poses = results[:10]

    with open(scores_file, 'w') as scores:
        scores.write("Molecule ID, First Docking, Second Docking\n")
        for second_score, pdbqt_file, best_pose_lines in top_10_poses:
            molecule_id = pdbqt_file[:-6]
            final_pose_path = os.path.join(final_poses_dir, pdbqt_file)

            # Write the best pose to a new file
            with open(final_pose_path, 'w') as file:
                file.writelines(best_pose_lines)

            # Convert to MOL2 using Open Babel
            mol2_file = final_pose_path.replace(".pdbqt", ".mol2")
            conversion_cmd = f"obabel {final_pose_path} -O {mol2_file}"
            conversion_result = subprocess.run(conversion_cmd, shell=True, stderr=subprocess.PIPE, text=True)

            # Check for conversion issues and update scores file
            if "Warning  in PerceiveBondOrders" in conversion_result.stderr:
                os.remove(final_pose_path)
                os.remove(mol2_file)
            else:
                os.remove(final_pose_path)
                first_docking_score, _ = get_best_score_and_pose(os.path.join(results_dir1, pdbqt_file))
                scores.write(f"{molecule_id}, {first_docking_score}, {second_score}\n")

    print("Analysis completed. Final poses and scores are saved in the specified directories.")

if __name__ == '__main__':
    # Example of how to call the function
    results_dir1 = "dock1_results_2bv6"
    results_dir2 = "dock2_results_2bv6"
    final_poses_dir = "final_poses_2bv6"
    scores_file = "scores_final_2bv6.txt"
    analyze_results(results_dir1, results_dir2, final_poses_dir, scores_file)
