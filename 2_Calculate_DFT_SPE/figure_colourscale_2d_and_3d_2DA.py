'''
This script is used to generate a 2d and a 3d plot from the conformer data 

Made by Houyang Xu and Paula Teeuwen, last updated by Houyang Xu 
Last updated: 12 Jul 2025
Adjusted for Github: 20 Nov 2025
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams

# Set font preferences
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica Neue', 'sans-serif']
rcParams['font.size'] = 6  # Base font size

# Read the data file
print("Reading conformer data...")
data = pd.read_csv('conformer_distances_data.csv')
print(f"Loaded {len(data)} conformers")

# Extract columns
dist_a = data['Avg_1-21_33-46']
dist_b = data['Avg_21-33_1-46']
energies = data['Relative_Energy_kcal']

# Filter data to only include energies up to 100 kcal/mol
mask = energies <= 100.0
dist_a_filtered = dist_a[mask]
dist_b_filtered = dist_b[mask]
energies_filtered = energies[mask]
print(f"Conformers with energy ≤ 100 kcal/mol: {len(energies_filtered)}")

# Create custom colormap with specified colors
dark_green = (32/255, 127/255, 76/255)
dark_red = (194/255, 31/255, 48/255)
colors = [dark_green, 'green', 'yellow', 'orange', dark_red]
n_bins = 100
cmap = mcolors.LinearSegmentedColormap.from_list('CustomGreenToRed', colors, N=n_bins)

# Normalize colors to 0-100 kcal/mol range
norm = mcolors.Normalize(vmin=0, vmax=100)

# ---------------------------------------------------------------------------
# 2D scatter plot
# ---------------------------------------------------------------------------
print("\nCreating 2D scatter plot...")

# Create figure (15 cm x 11 cm)
fig_width_2d = 15 / 2.54
fig_height_2d = 11 / 2.54
fig_2d, ax_2d = plt.subplots(figsize=(fig_width_2d, fig_height_2d))

# Calculate edge colors (0.7x RGB of face colors)
face_colors = cmap(norm(energies_filtered))
edge_colors = face_colors * 0.7
edge_colors[:, 3] = 1.0  # Keep alpha at 1

# Create 2D scatter plot
scatter_2d = ax_2d.scatter(dist_a_filtered, dist_b_filtered, 
                          c=energies_filtered, 
                          cmap=cmap,
                          norm=norm,
                          marker='o', 
                          s=15, 
                          alpha=0.8,
                          linewidth=0.5,
                          edgecolors=edge_colors)

# Set axis limits
ax_2d.set_xlim(14.0, 18.0)
ax_2d.set_ylim(6.25, 12.25)

# Set major and minor ticks
ax_2d.xaxis.set_major_locator(MultipleLocator(0.5))
ax_2d.xaxis.set_minor_locator(MultipleLocator(0.25))
ax_2d.yaxis.set_major_locator(MultipleLocator(1.0))
ax_2d.yaxis.set_minor_locator(MultipleLocator(0.5))

# Enable minor ticks and disable grid
ax_2d.minorticks_on()
ax_2d.grid(False)

# Set axis labels
ax_2d.set_xlabel('N···N distance $\it{l}$ (long edge, Å)', fontsize=6)
ax_2d.set_ylabel('N···N distance $\it{s}$ (short edge, Å)', fontsize=6)

# Add colorbar
cbar_2d = plt.colorbar(scatter_2d, ax=ax_2d)
cbar_2d.set_label('Relative Energy (kcal·mol⁻¹)', fontsize=6)
cbar_2d.set_ticks(np.arange(0, 101, 20))
cbar_2d.ax.tick_params(labelsize=6)

# Adjust tick label sizes
ax_2d.tick_params(axis='both', which='major', labelsize=6)
ax_2d.tick_params(axis='both', which='minor', labelsize=6)

# Tight layout and save
plt.tight_layout()
output_2d = 'conformer_scatter_plot.svg'
plt.savefig(output_2d, dpi=600, bbox_inches='tight', transparent=True)
print(f"2D plot saved as '{output_2d}'")

# ---------------------------------------------------------------------------
# 3D scatter plot (coloured version)
# ---------------------------------------------------------------------------
print("\nCreating 3D scatter plot (coloured)...")

# Create figure (18.39 cm x 13.79 cm)
fig_width_3d = 18.39 / 2.54
fig_height_3d = 13.79 / 2.54
fig_3d = plt.figure(figsize=(fig_width_3d, fig_height_3d))
ax_3d = fig_3d.add_subplot(111, projection='3d')

# Create 3D scatter plot
scatter_3d = ax_3d.scatter(dist_a_filtered, dist_b_filtered, energies_filtered,
                          c=energies_filtered, 
                          cmap=cmap,
                          norm=norm,
                          marker='o', 
                          s=20, 
                          alpha=0.7,
                          linewidth=0.3)

# Set axis limits
ax_3d.set_xlim(14.0, 18.0)
ax_3d.set_ylim(6.25, 12.25)
ax_3d.set_zlim(0, 100)

# Set major ticks
ax_3d.xaxis.set_major_locator(MultipleLocator(0.5))
ax_3d.yaxis.set_major_locator(MultipleLocator(1.0))
ax_3d.zaxis.set_major_locator(MultipleLocator(20))

# Set axis labels
ax_3d.set_xlabel('N···N distance $\it{l}$ (long edge, Å)', fontsize=9, labelpad=5)
ax_3d.set_ylabel('N···N distance $\it{s}$ (short edge, Å)', fontsize=9, labelpad=5)
ax_3d.set_zlabel('Relative Energy (kcal·mol⁻¹)', fontsize=9, labelpad=5)

# Add colorbar
cbar_3d = plt.colorbar(scatter_3d, ax=ax_3d, shrink=0.6, aspect=20)
cbar_3d.set_label('Relative Energy (kcal·mol⁻¹)', fontsize=9)
cbar_3d.set_ticks(np.arange(0, 101, 20))
cbar_3d.ax.tick_params(labelsize=9)

# Adjust tick label sizes and viewing angle
ax_3d.tick_params(axis='both', which='major', labelsize=9)
ax_3d.view_init(elev=20, azim=45)

# Tight layout and save
plt.tight_layout()
output_3d = 'conformer_3d_plot.svg'
plt.savefig(output_3d, dpi=600, bbox_inches='tight', transparent=True)
print(f"3D plot (coloured) saved as '{output_3d}'")

# ---------------------------------------------------------------------------
# 3D scatter plot (monochrome version)
# ---------------------------------------------------------------------------
# print("\nCreating 3D scatter plot (monochrome)...")
# 
# fig_3d_mono = plt.figure(figsize=(fig_width_3d, fig_height_3d))
# ax_3d_mono = fig_3d_mono.add_subplot(111, projection='3d')
# 
# # Single colour version
# scatter_3d_mono = ax_3d_mono.scatter(dist_a_filtered, dist_b_filtered, energies_filtered,
#                                     color='steelblue',
#                                     marker='o', 
#                                     s=15, 
#                                     alpha=0.6,
#                                     linewidth=0.3,
#                                     edgecolors='navy')
# 
# # Set axis limits
# ax_3d_mono.set_xlim(14.0, 18.0)
# ax_3d_mono.set_ylim(6.25, 12.25)
# ax_3d_mono.set_zlim(0, 100)
# 
# # Set major ticks
# ax_3d_mono.xaxis.set_major_locator(MultipleLocator(0.5))
# ax_3d_mono.yaxis.set_major_locator(MultipleLocator(1.0))
# ax_3d_mono.zaxis.set_major_locator(MultipleLocator(20))
# 
# # Set axis labels
# ax_3d_mono.set_xlabel('N···N distance $\it{l}$ (long edge, Å)', fontsize=9, labelpad=5)
# ax_3d_mono.set_ylabel('N···N distance $\it{s}$ (short edge, Å)', fontsize=9, labelpad=5)
# ax_3d_mono.set_zlabel('Relative Energy (kcal·mol⁻¹)', fontsize=9, labelpad=5)
# 
# # Adjust tick label sizes and viewing angle
# ax_3d_mono.tick_params(axis='both', which='major', labelsize=9)
# ax_3d_mono.view_init(elev=20, azim=45)
# 
# # Tight layout and save
# plt.tight_layout()
# output_3d_mono = 'conformer_3d_plot_mono.svg'
# plt.savefig(output_3d_mono, dpi=600, bbox_inches='tight', transparent=True)
# print(f"3D plot (monochrome) saved as '{output_3d_mono}'")

# ---------------------------------------------------------------------------
# Display statistics
# ---------------------------------------------------------------------------
print(f"\nPlot statistics:")
print(f"X-axis (distance $\it{l}$) range in data: {dist_a_filtered.min():.3f} - {dist_a_filtered.max():.3f} Å")
print(f"Y-axis (distance $\it{s}$) range in data: {dist_b_filtered.min():.3f} - {dist_b_filtered.max():.3f} Å")
print(f"Energy range in data: {energies_filtered.min():.3f} - {energies_filtered.max():.3f} kcal·mol⁻¹")

plt.show()