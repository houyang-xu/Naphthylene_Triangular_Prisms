'''
This script is used to generate a simplified data file

Made by Houyang Xu and Paula Teeuwen, last updated by Houyang Xu 
Last updated: 12 Jul 2025
Adjusted for Github: 20 Nov 2025
'''

import os
from ase.io import read
import numpy as np

# Load all structures from the trajectory file
print("Loading structures from combined_conformers_DFT_SPE.xyz...")
config_all = read('combined_conformers_DFT_SPE.xyz', index=":")
print(f"Total number of structures loaded: {len(config_all)}")

# Define atom indices for distance calculations (1-based indexing)
atom1, atom21, atom33, atom46 = 1, 21, 33, 46

# Conversion factor from eV to kcal/mol
ev_to_kcal = 23.060548867

# Find minimum energy first
min_energy = None
for structure in config_all:
    dft_energy = structure.info.get("DFT_energy")
    if dft_energy is not None:
        if min_energy is None or dft_energy < min_energy:
            min_energy = dft_energy

if min_energy is None:
    print("Error: No DFT energies found in file")
    exit(1)

print(f"Minimum energy found: {min_energy:.6f} eV")

# Open output file
with open("conformer_distances_data.csv", "w") as f:
    # Write header
    f.write("Sequence_Number,Conformer_Number,Avg_1-21_33-46,Avg_21-33_1-46,Relative_Energy_kcal\n")
    
    # Process each structure
    valid_count = 0
    for idx, structure in enumerate(config_all):
        positions = structure.get_positions()
        
        # Calculate individual distances (convert 1-based to 0-based indexing)
        dist_1_21 = np.linalg.norm(positions[atom1 - 1] - positions[atom21 - 1])
        dist_33_46 = np.linalg.norm(positions[atom33 - 1] - positions[atom46 - 1])
        dist_21_33 = np.linalg.norm(positions[atom21 - 1] - positions[atom33 - 1])
        dist_1_46 = np.linalg.norm(positions[atom1 - 1] - positions[atom46 - 1])
        
        # Calculate average distances
        avg_a = (dist_1_21 + dist_33_46) / 2.0  # Average of (1-21) and (33-46)
        avg_b = (dist_21_33 + dist_1_46) / 2.0  # Average of (21-33) and (1-46)
        
        # Extract conformer index from metadata
        conformer_idx = None
        for key in structure.info.keys():
            if key.endswith(':') and key[:-1].isdigit():
                conformer_idx = int(key[:-1])
                break
        
        # Extract DFT_energy and calculate relative energy
        dft_energy = structure.info.get("DFT_energy")
        if dft_energy is not None:
            relative_energy = (dft_energy - min_energy) * ev_to_kcal
            
            # Only include conformers with relative energy ≤ 100 kcal/mol
            if relative_energy <= 100.0:
                # Write data row
                sequence_num = idx + 1
                conformer_num = conformer_idx if conformer_idx is not None else sequence_num
                f.write(f"{sequence_num},{conformer_num},{avg_a:.6f},{avg_b:.6f},{relative_energy:.6f}\n")
                valid_count += 1
        else:
            print(f"Warning: Sequence {idx + 1}: DFT_energy not found, skipping.")

print(f"\nData file 'conformer_distances_data.csv' created successfully!")
print(f"Total conformers within 100 kcal/mol: {valid_count}")
print(f"Total conformers analyzed: {len(config_all)}")
print(f"Conformers excluded (>100 kcal/mol or no energy): {len(config_all) - valid_count}")