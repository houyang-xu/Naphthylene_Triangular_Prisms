'''
This script was used to create a 3D plot showing the DFT calculated single point energies versus the dihedral angles of the conformationally flexible quadrilateral linker. 
The z-axis is truncated at 200 kcal/mol. On the xy-plane a contourplot of the plot is created.

P.C.P. Teeuwen 
Last updated: 25-04-2025
Adjusted for GitHub: 19-11-2025
'''

from ase.io import read
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

# Load the multi-frame XYZ file
structures = list(read("combined_conformers_DFT_SPE.xyz", index=":"))

for i, atoms in enumerate(structures[:5]):
    print(f"Structure {i} info:")
    print(atoms.info)

data = []

for atoms in structures:
    info = atoms.info
    try:
        DA1 = float(info.get("DA1", "nan"))
        DA2 = float(info.get("DA2", "nan"))
        energy = float(info.get("DFT_energy", "nan"))
        if not any(map(np.isnan, [DA1, DA2, energy])):
            data.append([DA1, DA2, energy])
    except ValueError:
        continue

matrix = np.array(data)

min_energy = np.min(matrix[:, 2])
matrix[:, 2] -= min_energy
matrix[:, 2] *= 23.0609  # Convert energy units eV to kcal/mol

print("Matrix shape:", matrix.shape)
print("First few rows:\n", matrix[:5])

DA1, DA2, energy = matrix[:, 0], matrix[:, 1], matrix[:, 2]
z_cutoff = 200  # Set z-axis cutoff to 200 kcal/mol

grid_res = 100
xi = np.linspace(min(DA1), max(DA1), grid_res)
yi = np.linspace(min(DA2), max(DA2), grid_res)
Xi, Yi = np.meshgrid(xi, yi)
    
Zi = griddata((DA1, DA2), energy, (Xi, Yi), method='cubic')

# Handle NaNs by setting them to a fixed value or masking them
Zi_masked = np.ma.array(Zi, mask=np.isnan(Zi))
Zi_clipped = np.where(Zi_masked > z_cutoff, z_cutoff, Zi_masked)  # Limit Z values to 200

# Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(Xi, Yi, Zi_clipped, cmap='RdYlGn_r', edgecolor='none', alpha=0.75)
ax.contour(Xi, Yi, Zi_clipped, zdir='z', cmap='RdYlGn_r', alpha=0.6, offset=-50)  # Add contour plot at z=0

ax.set_xlabel('DA1 (°)')
ax.set_ylabel('DA2 (°)')
ax.set_zlabel('DFT Energy (kcal/mol)')
ax.set_title('DFT SPE Energy Surface')

# Set z-limits to 0-200
ax.set_zlim(-50, 200)
ax.grid(False) 

fig.colorbar(surf, ax=ax, label='DFT Energy (kcal/mol)')  # Colorbar for DFT Energy values
plt.tight_layout()
plt.savefig("DA1_DA2_DFT_surface_contour_0_200.png", dpi=300)
plt.show()
