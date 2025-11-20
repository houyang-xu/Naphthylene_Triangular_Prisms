
'''
This script was used to append all atoms-dft.xyz files created by ASE after the DFT single-point energy calculations.
The ORCA output files were stored in scratch, however for the Github they are placed in the folder "Orca".
The new file with all the combined conformers is saved as "combined_conformers_DFT_SPE.xyz"

P.C.P. Teeuwen
Last updated: 11-12-2024
Adjusted for GitHub: 19-11-2025
'''

import os
import argparse

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Append DFT SPE results")
parser.add_argument('--DFTfolder', type=str, help='Folder containing the results for the SPE DFT calculations on all the conformers for a certain input structure')
args = parser.parse_args()
DFT_folder = args.DFTfolder
os.chdir(DFT_folder)

# Output file name
output_filename = f"combined_conformers_DFT_SPE.xyz"

# Open the output file for writing
with open(output_filename, 'w') as outfile:
    # List all directories in the DFT_results directory
    orca_job_folders = [d for d in os.listdir(DFT_folder) if d.startswith("ORCA_job_") and os.path.isdir(os.path.join(DFT_folder, d))]
    
    # Sort the directories by the numeric suffix
    orca_job_folders.sort(key=lambda x: int(x.split("_")[-1]))
    
    first_file = True
    
    # Loop through each ORCA_job_i folder in order
    for folder in orca_job_folders:
        filepath = os.path.join(DFT_folder, folder, "atoms-dft.xyz")
        if os.path.exists(filepath):
            # Read the contents of the file
            with open(filepath, 'r') as infile:
                content = infile.read().strip()
                content = content.replace("mp2_forces", "DFT_forces")
                content = content.replace("mp2_energy", "DFT_energy")                
                # Write content to the output file
                if not first_file:
                    outfile.write('\n')  # Separate files with a newline, not just blank lines
                outfile.write(content)
                first_file = False

print(f"All atoms-dft.xyz files have been combined into {output_filename}")

