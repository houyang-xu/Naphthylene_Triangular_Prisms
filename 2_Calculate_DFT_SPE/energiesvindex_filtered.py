'''
P.C.P. Teeuwen
Last updated: 11-12-2024
Adjusted for GitHub: 19-11-2025
'''

import os
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib import colormaps

# Lists to store distances for each pair of atoms
config_all = read('combined_conformers_DFT_SPE.xyz', index=":")
energies = []

# total_indices = 1296

#Remove doubles by going from 334-1000
# Conformer 334: DA1=0, DA2=0 
# Conformer 1000: DA1=180, DA2=180
start_conformer = 334
end_conformer = 1000
start_index = start_conformer - 1
end_index = end_conformer - 1


# Cannot simplify by taking energies directly from the combined conformers since we want to find those calculations that failed and identify them in the plot.
for idx in range(start_index, end_index + 1):
    file_path = f"ORCA/ORCA_job_{idx}/atoms-dft.xyz"
    
    if os.path.exists(file_path):
        # Read the file if it exists
        structure = read(file_path)
        positions = structure.get_positions()

        # Extract DFT_energy from structure metadata
        dft_energy = structure.info.get("DFT_energy")
        if dft_energy is not None:
            energies.append(dft_energy)
        else:
            energies.append(None)
            print(f"Index {idx}: DFT_energy not found in metadata.")
    else:
        # File does not exist
        energies.append(None)
        print(f"Index {idx}: File '{file_path}' does not exist.")

# Filter out None values to locate the missing indices
indices = np.arange(start_conformer, end_conformer + 1)  # Use 1-based numbering for clarity
valid_indices = [i for i, energy in enumerate(energies) if energy is not None]
valid_energies = [energy for energy in energies if energy is not None]
missing_indices = [i for i, energy in enumerate(energies) if energy is None]

print(' ')
print('missing indices are:')
print(missing_indices)

min_energy = np.min(valid_energies)
shifted_valid_energies = (valid_energies - min_energy)*23.0609 #Convert energy units from eV to kcal/mol
energy_threshold = 200

# Create the plot, including vertical lines to indicate missing values
plt.figure(figsize=(10, 6))
plt.plot(valid_indices, shifted_valid_energies, linestyle="-", color="blue", alpha=0.7)

label_added = False  # flag to track if label is already added

for idx in missing_indices:
    if not label_added:
        plt.axvline(x=idx, color="grey", linestyle="--", linewidth=1.5, alpha=0.8,
                    label="No DFT energy calculated", zorder=1)
        label_added = True
    else:
        plt.axvline(x=idx, color="grey", linestyle="--", linewidth=1.5, alpha=0.8, zorder=1)

plt.axhline(y=energy_threshold, color="orange", linestyle="--", linewidth=2, label=f"Threshold: {energy_threshold:.2f}", zorder=2)

# Add labels and title
plt.xlabel("Index of Conformer")
plt.ylabel("DFT Energy (kcal/mol)")
plt.title("DFT Energy vs. Structure Index")

# Highlight the minimum and maximum energies
min_energy = np.min(shifted_valid_energies) #should be 0
max_energy = np.max(shifted_valid_energies)
min_index = valid_indices[np.argmin(shifted_valid_energies)]
max_index = valid_indices[np.argmax(shifted_valid_energies)]

plt.scatter(min_index, min_energy, color="green", label=f"Min Energy: {min_energy:.2f}", zorder=5)
plt.scatter(max_index, max_energy, color="red", label=f"Max Energy: {max_energy:.2f}", zorder=5)

plt.legend()

# Save the figure as a PNG file
plt.savefig("energiesvindex.png", format="png", dpi=300)  # dpi=300 for high quality




