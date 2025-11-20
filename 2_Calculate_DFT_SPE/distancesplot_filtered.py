'''
P.C.P. Teeuwen
Last updated: 21-02-2024
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

# Load all structures from the trajectory file
config_all = read('combined_conformers_DFT_SPE.xyz', index=":")
#doesn't include the conformers it couldn't find the energy for


# Provide the indices of the N-N distances that should be recorded. Use a 1-based indexing here.
# Pymol uses 1-based indexing, which is converted into 0-based indexing as required by ASE 
atoma1, atoma2 = 1, 21
atomb1, atomb2 = 21, 33

nitrogen1 = 1
nitrogen2 = 21
nitrogen3 = 33
nitrogen4 = 46

# Lists to store energies of the conformers and distances for each pair of atoms
energies = []
distancesa = [] 
distancesb = []
#distancesc = []

# Access the first structure to check what atoms are located from the given indexes
first_structure = config_all[0]  
print(f"Atom pair A: {first_structure[nitrogen1 - 1].symbol} (index {nitrogen1}) - {first_structure[nitrogen2 - 1].symbol} (index {nitrogen2})")
print(f"Atom pair B: {first_structure[nitrogen2 - 1].symbol} (index {nitrogen2}) - {first_structure[nitrogen3 - 1].symbol} (index {nitrogen3})")
print(f"Atom pair C: {first_structure[nitrogen3 - 1].symbol} (index {nitrogen3}) - {first_structure[nitrogen4 - 1].symbol} (index {nitrogen4})")#
print(f"Atom pair D: {first_structure[nitrogen4 - 1].symbol} (index {nitrogen4}) - {first_structure[nitrogen1 - 1].symbol} (index {nitrogen1})")#

# Here the doubles are removed and only the results for a selection of conformers is shown
start_conformer = 334 
end_conformer = 1000
#start_conformer = 0
#end_conformer = len(config_all) 

start_index = start_conformer - 1
end_index = end_conformer - 1

# Loop through each structure to extract distances and energies
for idx in range(start_index, end_index + 1):
    structure = config_all[idx]
    positions = structure.get_positions()
    # Calculate distances for atom pairs
    dista = (np.linalg.norm(positions[nitrogen1 - 1] - positions[nitrogen2 - 1]) + np.linalg.norm(positions[nitrogen3 - 1] - positions[nitrogen4 - 1]))/2
    distancesa.append(dista)

    distb = (np.linalg.norm(positions[nitrogen2 - 1] - positions[nitrogen3 - 1]) + np.linalg.norm(positions[nitrogen4 - 1] - positions[nitrogen1 - 1]))/2
    distancesb.append(distb)

    # Extract DFT_energy from structure metadata
    dft_energy = structure.info.get("DFT_energy")
    if dft_energy is not None:
        energies.append(dft_energy)
    else:
            #If the SPE calculation failed a None is added, so that those indices are still plotted / indices still correspond correctly
        energies.append(None)
        print(f"Index {idx}: DFT_energy not found in metadata.")

print(f"Number of structures: {len(config_all)}")

# Convert lists to numpy arrays for easier manipulation
distancesa = np.array(distancesa)
distancesb = np.array(distancesb)
energies = np.array(energies)

# From the list with energies, find the minimum energy and shift everything such that the lowest 
# value is 0. Additionally convert eV values (DFT_energy in combined xyz file are from ASE's atoms-dft.xyz and are given in eV) into kcal/mol.

energies_nonone = [e for e in energies if e is not None]
min_energy = np.nanmin(energies_nonone)
energy_threshold = 100.01
shifted_energies = np.array(
    [(e - min_energy) * 23.0609 if e is not None else None for e in energies]) #conversion eV to kcal/mol

# Apply a "mask" where only conformers with an energy below the threshold are plotted
mask = shifted_energies < energy_threshold
filtered_distancesa = distancesa[mask]
filtered_distancesb = distancesb[mask]
filtered_energies = shifted_energies[mask]

# Normalize energy values for colormap
norm = Normalize(vmin=filtered_energies.min(), vmax=filtered_energies.max())
cmap = colormaps["RdYlGn_r"]

min_a, max_a = filtered_distancesa.min(), filtered_distancesa.max()
min_b, max_b = filtered_distancesb.min(), filtered_distancesb.max()

print(f"Atom pair A ({atoma1}-{atoma2}): Min distance = {min_a:.2f} Å, Max distance = {max_a:.2f} Å")
print(f"Atom pair B ({atomb1}-{atomb2}): Min distance = {min_b:.2f} Å, Max distance = {max_b:.2f} Å")

# Plot the correlation between distancesa and distancesb
scatter = plt.scatter(
    filtered_distancesa, filtered_distancesb, c=filtered_energies, cmap="RdYlGn_r", alpha=0.7, edgecolor="k")
plt.colorbar(scatter, label="DFT Energy (kcal/mol)")  # Add a colorbar for reference

# Allows adding the rectangular boxes in the next step
ax = plt.gca()  # Get the current axis

plt.xlabel(f"N-N distance a+c/2 (Å)")
plt.ylabel(f"N-N distance b+d/2 (Å)")
plt.title("Correlation between short (a,c) and long (b,d) N-N distances")
plt.grid(True)

# Save the figure as a PNG file
plt.savefig("distances_plot_filtered.png", format="png", dpi=300)  # dpi=300 for high quality

# Display the plot
plt.show()



