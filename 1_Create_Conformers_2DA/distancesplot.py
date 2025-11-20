'''
This scirpt is used to extract the N-N distances from the generated conformer database and plot the long distance vs the short distance to visualise the range.

P.C.P. Teeuwen
Last updated: 11-12-2024
Adjusted for Github: 19-11-2025
'''

import os
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Load all structures from the trajectory file
config_all = read('Conformers_2DA.xyz', index=":")

# Provide 1-based indexing here (e.g. from pymol)
# The code adjusts this later to 0-based indexing
atoma1, atoma2 = 1, 21
atomb1, atomb2 = 21, 33

# Lists to store distances for each pair of atoms
distancesa = [] 
distancesb = []
#distancesc = []

first_structure = config_all[0]  # Access the first structure
print(f"Atom pair A: {first_structure[atoma1 - 1].symbol} (index {atoma1}) - {first_structure[atoma2 - 1].symbol} (index {atoma2})")
print(f"Atom pair B: {first_structure[atomb1 - 1].symbol} (index {atomb1}) - {first_structure[atomb2 - 1].symbol} (index {atomb2})")
#print(f"Atom pair C: {first_structure[atomc1 - 1].symbol} (index {atomc1}) - {first_structure[atomc2 - 1].symbol} (index {atomc2})")#

# Calculate distances for each structure in the trajectory
for structure in config_all:
    positions = structure.get_positions()
    
    dista = np.linalg.norm(positions[atoma1 - 1 ] - positions[atoma2 - 1])
    distancesa.append(dista)

    distb = np.linalg.norm(positions[atomb1 - 1] - positions[atomb2 - 1])
    distancesb.append(distb)

    #distc = np.linalg.norm(positions[atomc1 - 1] - positions[atomc2 - 1])
    #distancesc.append(distc)

# Convert lists to numpy arrays (optional, for easier manipulation)
distancesa = np.array(distancesa)
distancesb = np.array(distancesb)
#distancesc = np.array(distancesc)

min_a, max_a = distancesa.min(), distancesa.max()
min_b, max_b = distancesb.min(), distancesb.max()
#min_c, max_c = distancesc.min(), distancesc.max()

print(f"Atom pair A ({atoma1}-{atoma2}): Min distance = {min_a:.2f} Å, Max distance = {max_a:.2f} Å")
print(f"Atom pair B ({atomb1}-{atomb2}): Min distance = {min_b:.2f} Å, Max distance = {max_b:.2f} Å")
#print(f"Atom pair C ({atomc1}-{atomc2}): Min distance = {min_c:.2f} Å, Max distance = {max_c:.2f} Å")

# Plot the correlation between distancesa and distancesb
plt.figure(figsize=(8, 6))
plt.scatter(distancesa, distancesb, alpha=0.5, color="b")

# Adding rectangular boxes
ax = plt.gca()  # Get the current axis

# Green box coordinates with max and min values for experimental result 2
green_box_x = (15.186, 15.536)  # x range
green_box_y = (10.300, 10.929)  # y range
green_box_width = green_box_x[1] - green_box_x[0]
green_box_height = green_box_y[1] - green_box_y[0]

green_box = Rectangle(
    (green_box_x[0], green_box_y[0]), 
    green_box_width, 
    green_box_height,
    linewidth=1.5,
    edgecolor='green',
    facecolor='none'
)
ax.add_patch(green_box)

# Red box coordinates with max and min values for experimental result 1
red_box_x = (15.410, 16.175)  # x range
red_box_y = (9.095, 10.360)   # y range
red_box_width = red_box_x[1] - red_box_x[0]
red_box_height = red_box_y[1] - red_box_y[0]

red_box = Rectangle(
    (red_box_x[0], red_box_y[0]), 
    red_box_width, 
    red_box_height,
    linewidth=1.5,
    edgecolor='red',
    facecolor='none'
)
ax.add_patch(red_box)

# Box with average distances for trigonal prism 1
largestTPa = [15.926]
largestTPb = [9.698]
plt.scatter(largestTPa, largestTPb, alpha=1, color="r")

# Box with average distances for trigonal prism 2
smallestTPa = [15.331]
smallestTPb = [10.625]
plt.scatter(smallestTPa, smallestTPb, alpha=1, color="g")

plt.xlabel(f"N-N distance a (Å)")
plt.ylabel(f"N-N distance b (Å)")
plt.title("Correlation between N-N distances a and b")
plt.grid(True)

# Save the figure as a PNG file
plt.savefig("distances_plot.png", format="png", dpi=300)  # dpi=300 for high quality

# Display the plot
plt.show()


