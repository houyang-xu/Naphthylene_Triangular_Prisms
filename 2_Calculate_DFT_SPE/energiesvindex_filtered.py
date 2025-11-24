'''
<<<<<<< HEAD
This script is used to generate a mp4 video or gif from the conformer data 

Made by Houyang Xu and Paula Teeuwen, last updated by Houyang Xu 
Last updated: 12 Jul 2025
Adjusted for Github: 21 Nov 2025
=======
P.C.P. Teeuwen
Last updated: 11-12-2024
Adjusted for GitHub: 19-11-2025
>>>>>>> 816dce403ab5b672a91e6d374518c5f847eec17f
'''

import os
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib import colormaps

<<<<<<< HEAD
# Set figure properties
plt.rcParams.update({
    'font.size': 9,                      # Set all font sizes to 9
    'font.family': 'sans-serif',                              
    'font.sans-serif': ['Helvetica Neue', 'Arial'],  
    'axes.labelsize': 9,                  # Axis labels
    'xtick.labelsize': 9,                  # Tick labels
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
})

# Lists to store distances for each pair of atoms
config_all = read('Hcombined_conformers_DFT_SPE.xyz', index=":")
energies = []

total_indices = 1296
=======
# Lists to store distances for each pair of atoms
config_all = read('combined_conformers_DFT_SPE.xyz', index=":")
energies = []

# total_indices = 1296
>>>>>>> 816dce403ab5b672a91e6d374518c5f847eec17f

#Remove doubles by going from 334-1000
# Conformer 334: DA1=0, DA2=0 
# Conformer 1000: DA1=180, DA2=180
<<<<<<< HEAD
# Here the doubles are removed and only the results for a selection of conformers is shown
#start_conformer = 334
#end_conformer = 1000
start_conformer = 0
end_conformer = len(config_all)

start_index = start_conformer - 1
end_index = end_conformer - 1

for idx in range(start_index, end_index + 1):
    file_path = f"ORCA/ORCA_job_{idx}/atoms-dft.xyz"

=======
start_conformer = 334
end_conformer = 1000
start_index = start_conformer - 1
end_index = end_conformer - 1


# Cannot simplify by taking energies directly from the combined conformers since we want to find those calculations that failed and identify them in the plot.
for idx in range(start_index, end_index + 1):
    file_path = f"ORCA/ORCA_job_{idx}/atoms-dft.xyz"
    
>>>>>>> 816dce403ab5b672a91e6d374518c5f847eec17f
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

<<<<<<< HEAD

# `energies` is a list that may contain None values
indices = np.arange(start_conformer, end_conformer + 1)  # Use 1-based numbering for clarity

# Filter out None values to locate the missing indices
valid_indices = [i for i, energy in enumerate(energies) if energy is not None]
valid_energies = [energy for energy in energies if energy is not None]
missing_indices = [i for i, energy in enumerate(energies) if energy is None]
print(missing_indices)

min_energy = np.min(valid_energies)
index_min_energy = np.argmin(valid_energies)
print(f'Lowest conformer index is: {index_min_energy}')
shifted_valid_energies = (valid_energies - min_energy)*23.060

print(f"Min shifted energy: {np.min(shifted_valid_energies)}")
print(f"Max shifted energy: {np.max(shifted_valid_energies)}")

# energy_threshold = 1000
energy_threshold = 100

# Create the plot, including vertical lines to indicate missing values
plt.figure(figsize=(7.2, 4.05))
plt.plot(valid_indices, shifted_valid_energies, linestyle="-", color="blue", alpha=0.7)

for idx in missing_indices:
    plt.axvline(x=idx, color="grey", linestyle="--", linewidth=1, alpha=0.6, zorder=1)

plt.axhline(y=energy_threshold, color="orange", linestyle="--", linewidth=1, 
           label=f"Threshold: {energy_threshold:.2f} kcal·mol⁻¹", zorder=2, alpha=0.6)


# Add labels and title
plt.xlabel("Index of Structure")
plt.ylabel("DFT Energy (kcal·mol⁻¹)")
plt.title("DFT Energy vs. Structure Index")
# Add this line before saving the figure
plt.ylim(-10, 200)  # Sets y-axis limits from 0 to 110 kcal/mol
=======
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
>>>>>>> 816dce403ab5b672a91e6d374518c5f847eec17f

# Highlight the minimum and maximum energies
min_energy = np.min(shifted_valid_energies) #should be 0
max_energy = np.max(shifted_valid_energies)
min_index = valid_indices[np.argmin(shifted_valid_energies)]
max_index = valid_indices[np.argmax(shifted_valid_energies)]

<<<<<<< HEAD
plt.scatter(min_index, min_energy, color="green",s=10, label=f"Min Energy: {min_energy:.2f}", zorder=5)
# plt.scatter(max_index, max_energy, color="red",s=10, label=f"Max Energy: {max_energy:.2f}", zorder=5)

# Add a legend
plt.legend()

# Save the figure as a SVG file
plt.savefig("energiesvindex.svg", format="svg", dpi=400, transparent=True)  # dpi=00 for high quality
=======
plt.scatter(min_index, min_energy, color="green", label=f"Min Energy: {min_energy:.2f}", zorder=5)
plt.scatter(max_index, max_energy, color="red", label=f"Max Energy: {max_energy:.2f}", zorder=5)

plt.legend()

# Save the figure as a PNG file
plt.savefig("energiesvindex.png", format="png", dpi=300)  # dpi=300 for high quality

>>>>>>> 816dce403ab5b672a91e6d374518c5f847eec17f



