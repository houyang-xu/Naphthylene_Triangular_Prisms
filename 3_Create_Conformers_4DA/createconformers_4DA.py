'''
With this script a database is created for the conformers of a quadrilateral linker by iterating over four dihedral angles.
A checkfile is created "distances_5.csv" which prevents rerunning the screening after creation.
The file continues by plotting three graphs showing the range of N-N distances the linker can achieve as well as which parts of that range were accessed experimentally.

P.C.P. Teeuwen
Lasted updated by Houyang Xu
Last updated: 1 Aug 2025
Updated for GitHub: 24 Nov 2025
'''

from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
import numpy as np
import pickle
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.text as mtext
import matplotlib.font_manager as fm
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
import matplotlib.colors as mcolors
from matplotlib import rcParams

# Set font preferences with fallback
rcParams['font.family'] = 'sans-serif'
# rcParams['font.sans-serif'] = ['Helvetica Neue', 'Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
rcParams['font.sans-serif'] = ['Helvetica Neue', 'Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
rcParams['font.size'] = 6

# Create font properties for bold and bold italic
regular_font = fm.FontProperties(family='sans-serif', weight='normal')
bold_font = fm.FontProperties(family='sans-serif', weight='bold')
italic_font = fm.FontProperties(family='sans-serif', style='italic')
# bold_italic_font = fm.FontProperties(family='Helvetica', weight='bold', style='oblique')

# Use exact font files in MathText
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = regular_font.get_name()
mpl.rcParams['mathtext.bf'] = bold_font.get_name()  
mpl.rcParams['mathtext.default'] = 'bf'  
mpl.rcParams['text.usetex'] = False  
mpl.rcParams['font.family'] = regular_font.get_name()  
mpl.rcParams['font.size'] = 6  # Small font size

# ASE uses 0-based indexing
# Pymol uses 1-based indexing

nitrogen1 = 1
nitrogen2 = 21
nitrogen3 = 33
nitrogen4 = 46

def smiles_to_xyz_with_locked_dihedrals(smiles: str, output_filename: str, dihedral_dict: dict):
    """
    Generates a single XYZ file containing conformations with locked dihedral angles.

    Parameters:
        smiles (str): SMILES string of the molecule.
        output_filename (str): Name of the output XYZ file to save all conformers.
        dihedral_dict (dict): Dictionary specifying dihedral angles and ranges in the format:
            {
                "DA1": [index_A, index_B, index_C, index_D, [angle_1, angle_2, ..., angle_n]],
                ...
            }
            Note: Atom indices should be 1-based in input as obtained from e.g. Pymol and will be adjusted to 0-based.
    """

    distancesa = []
    distancesb = []
    distancesc = []
    distancesd = []
    # Generate the base molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol, addCoords=True)  # Adds protons to molecule
    AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Creates 3D molecule

    # Extract ranges for independent dihedrals (DA1/DA4 and DA2/DA3 are locked to create a C2 symmetric molecule)
#    da1_da4_angles = dihedral_dict["DA1"][4]  # Shared range for DA1 and DA4
#    da2_da3_angles = dihedral_dict["DA2"][4]  # Shared range for DA2 and DA3
    da1_angles = dihedral_dict["DA1"][4]
    da2_angles = dihedral_dict["DA2"][4]
    da3_angles = dihedral_dict["DA3"][4]
    da4_angles = dihedral_dict["DA4"][4]

    angle_combinations = itertools.product(da1_angles, da2_angles, da3_angles, da4_angles) # Create all possible pairs of angles
    print(len(da1_angles)*len(da2_angles)*len(da3_angles)*len(da4_angles))
    #return
    # Open the output file for writing all conformers
    with open(output_filename, 'w') as xyz_file:
        for conf_id, (angle1, angle2, angle3, angle4) in enumerate(angle_combinations):
            conf = mol.GetConformer()
            # Set DA1 to angle1
            for key in ["DA1"]:
                atom_indices = [idx - 1 for idx in dihedral_dict[key][:4]]  # Adjust to 0-based
                AllChem.SetDihedralDeg(conf, *atom_indices, float(angle1))

            for key in ["DA2"]:
                atom_indices = [idx - 1 for idx in dihedral_dict[key][:4]]  # Adjust to 0-based
                AllChem.SetDihedralDeg(conf, *atom_indices, float(angle2))

            for key in ["DA3"]:
                atom_indices = [idx - 1 for idx in dihedral_dict[key][:4]]  # Adjust to 0-based
                AllChem.SetDihedralDeg(conf, *atom_indices, float(angle3))

            for key in ["DA4"]:
                atom_indices = [idx - 1 for idx in dihedral_dict[key][:4]]  # Adjust to 0-based
                AllChem.SetDihedralDeg(conf, *atom_indices, float(angle4))

            # Write down four N···N distances. This time making them from an rdkit molecule instead of an ase object like in the other scripts. To
            # circumvent having to write a huge combined xyz file

            pos_n1 = np.array(conf.GetAtomPosition(nitrogen1 - 1))
            pos_n2 = np.array(conf.GetAtomPosition(nitrogen2 - 1))
            pos_n3 = np.array(conf.GetAtomPosition(nitrogen3 - 1))
            pos_n4 = np.array(conf.GetAtomPosition(nitrogen4 - 1))

            dista = np.linalg.norm(pos_n1 - pos_n2)
            distancesa.append(dista)

            distb = np.linalg.norm(pos_n2 - pos_n3)
            distancesb.append(distb)

            distc = np.linalg.norm(pos_n3 - pos_n4)
            distancesc.append(distc)

            distd = np.linalg.norm(pos_n4 - pos_n1)
            distancesd.append(distd)


            # Write the resulting conformation to the XYZ file, don't need to do that, just get distances
            #xyz_block = Chem.MolToXYZBlock(mol)
            #atoms = xyz_block.splitlines()[2:]

            #xyz_file.write(f"{len(atoms)}\n")  # Number of atoms
            #xyz_file.write(f" Conformer {conf_id + 1}: DA1={angle1}, DA2={angle2}, DA3={angle3}, DA4={angle4}\n")  # Comment line
            #xyz_file.write("\n".join(atoms) + "\n")

            #print(f"Saved conformer {conf_id + 1} to {output_filename}")
            #print(f"Saved distances belonging to {conf_id + 1}")

        return distancesa, distancesb, distancesc, distancesd

# Example usage
angles_dict = {
    "DA1": [38, 39, 40, 50, range(0, 360, 5)],  # Atom indices for dihedral 1
    "DA2": [13, 12, 11, 10, range(0, 360, 5)],  # Atom indices for dihedral 2
    "DA3": [38, 26, 27, 28, range(0, 360, 5)],  # Atom indices for dihedral 3
    "DA4": [13, 14, 15, 25, range(0, 360, 5)],  # Atom indices for dihedral 4
}

smiles = "NC1=CC(C=C2)=C(C=C1)C=C2C3=CC(C4=CC=C(C=C(N)C=C5)C5=C4)=C(C6=CC(C=CC(N)=C7)=C7C=C6)C=C3C8=CC=C(C=C(N)C=C9)C9=C8"

def load_and_sample_csv(filename, sample_frac=0.01):
    """Load large CSV in chunks and sample"""
    chunks = []
    for chunk in pd.read_csv(filename, chunksize=10000):  # 10000 rows at a time
        sampled_chunk = chunk.sample(frac=sample_frac)
        chunks.append(sampled_chunk)
    
    df = pd.concat(chunks, ignore_index=True)
    return df

if os.path.exists('distances_5.csv'):
#      df = load_and_sample_csv('distances_5.csv', sample_frac=0.1)  # Down sampling for speed. this keeps 10% = ~1040k points. This basically doesnt look different to using all points. Good when a quick view is needed
  df = load_and_sample_csv('distances_5.csv', sample_frac=1)  # Keep all points. Do this for final graphics.
    # Extract the columns back into individual arrays
    distancesa = df['DistancesA'].values
    distancesb = df['DistancesB'].values
    distancesc = df['DistancesC'].values
    distancesd = df['DistancesD'].values

    print("Data loaded successfully.")
else:

    distancesa, distancesb, distancesc, distancesd = smiles_to_xyz_with_locked_dihedrals(smiles, "HXLinker1_fourangles_all_conformers.xyz", angles_dict)
    print(f'Length of distances: {len(distancesa)}')
    print(f'Length of distances: {len(distancesb)}')
    print(f'Length of distances: {len(distancesc)}')
    print(f'Length of distances: {len(distancesd)}')

    # Combine your distance arrays into a single 2D array
    data = np.column_stack([distancesa, distancesb, distancesc, distancesd])

    # Convert to DataFrame and save as CSV
    df = pd.DataFrame(data, columns=['DistancesA', 'DistancesB', 'DistancesC', 'DistancesD'])
    df.to_csv('distances_5.csv', index=False)

alldistances = np.concatenate([distancesa, distancesb, distancesc, distancesd])

data = [distancesa, distancesb, distancesc, distancesd]

# --- Modified Violin Plot ---

# Create a figure for the violin plot
fig_violin, ax_violin = plt.subplots(figsize=(3.9, 2.7))

# Reduce the gap between violins by placing them closer together
# Original positions were [1, 2, 3, 4], now using smaller intervals
positions = [1, 1.7, 2.4, 3.1]  # Reduced spacing between positions

# Create the violin plot with reduced width to prevent overlap
violin_parts = ax_violin.violinplot(dataset=data,
                                    positions=positions,
                                    showmedians=True,
                                    widths=0.5)  # Reduced width from 0.75 to 0.5

# Set axis labels for the violin plot
ax_violin.set_xlabel('Distance type')
ax_violin.set_ylabel('N···N distance (Å)')

# Set x-axis ticks and labels with italic formatting
ax_violin.set_xticks(positions)
# Use LaTeX-style formatting for italics in matplotlib
ax_violin.set_xticklabels(['a', 'b', 'c', 'd'], fontproperties=italic_font, size=6)

# Adjust x-axis limits to center the plot nicely
ax_violin.set_xlim(0.5, 3.6)

# Calculate the median values for each dataset and add them as labels
medians = [np.median(distancesa), np.median(distancesb), np.median(distancesc), np.median(distancesd)]
for i, (pos, median) in enumerate(zip(positions, medians)):
    ax_violin.text(pos, median + 0.05, f'Median: {median:.2f}', 
                   horizontalalignment='center', 
                   verticalalignment='bottom', 
                   color='black',
                   fontsize=6)  # Slightly smaller font for median labels

# Optional: Add vertical grid lines for better readability
ax_violin.yaxis.grid(True, alpha=0.3)

# Save the violin plot as a SVG file
fig_violin.savefig('violin_plot_5.svg', dpi=600, bbox_inches='tight', transparent=True)

# ----- Violin Plot with datapoints -----------------

#smallestTP
data_A_1 = [15.4514,10.9164, 15.1883, 10.2954]
data_A_2 = [15.5376, 10.7482, 15.2594, 10.4751]
data_A_3 = [15.3112, 10.9205, 15.2528, 10.3914]
data_A = data_A_1 + data_A_2 + data_A_3

#midTP
data_B_1 = [16.0221, 9.6617, 15.8152, 9.4057]
data_B_2 = [14.4423, 12.1994, 14.4077, 10.6144]
data_B_3 = [14.886, 12.0502, 14.4883, 10.255]
data_B = data_B_1 + data_B_2 + data_B_3

#largestTP
data_C_a_1 = [16.1086,	9.5407,	15.9601, 9.2241]
data_C_a_2 = [15.8395,	10.0912, 15.6319, 9.832]
data_C_a_3 = [15.8461,	11.072,	15.3791, 9.4627]
data_C_b_1 = [15.9106,	9.9405,	15.7225, 9.6694]
data_C_b_2 = [16.0486,	9.3714,	15.9772, 9.3609]
data_C_b_3 = [16.1358,	9.8791,	15.9161, 9.2164]
data_C_c_1 = [16.2014,	9.4222,	16.1856, 9.1831]
data_C_c_2 = [16.0851,	10.3633, 15.6642, 9.868]
data_C_c_3 = [16.1063,	10.0081, 15.8883, 9.0993]
data_C_a = data_C_a_1 + data_C_a_2 + data_C_a_3
data_C_b = data_C_b_1 + data_C_b_2 + data_C_b_3
data_C_c = data_C_c_1 + data_C_c_2 + data_C_c_3
data_C = data_C_a + data_C_b + data_C_c

def density_based_jitter(energies, scale=0.1):
    # Calculate the kernel density estimate (KDE) for the energies
    kde = gaussian_kde(energies)
    density = kde.evaluate(energies)  # This returns the density at each energy point

    # The jitter should be inversely proportional to the density
    # High density -> high jitter; low density -> low jitter
    jitter_value = np.interp(density, (density.min(), density.max()), (0.0, scale))

    # Now we add a small random normal distribution to create the jitter
    category_jitter = np.random.normal(0, jitter_value, len(energies))
    return category_jitter

jitter_A = density_based_jitter(data_A, scale=0.1)
jitter_B = density_based_jitter(data_B, scale=0.1)
jitter_C = density_based_jitter(data_C, scale=0.1)

color_A = (0.725, 0.082, 0.227)  # Red tone
color_B = (0.0667, 0.4667, 0.6902)  # Blue tone
color_C = (0, 0.357, 0.188)  # Green tone

# Create figure and axes
fig_violin2, ax_violin2 = plt.subplots(figsize=(8.1768, 4.5995))

# Define x positions for data points
x_positions = [0, 1, 2, 3]  # Corresponding to the four violins

# Create violin plots
violin_parts2 = ax_violin2.violinplot([alldistances, alldistances, alldistances, alldistances], positions=[0,1,2,3], showmedians=False)

# Set different colors for each violin
colors = ['grey', color_A, color_B, color_C]  # Specify the colors

# Access each violin part and set its color
for i, pc in enumerate(violin_parts2['bodies']):
    pc.set_facecolor(colors[i])  # Set the color of the violin body
    pc.set_edgecolor(colors[i])    # Optional: Set the edge color
    pc.set_alpha(0.5)            # Optional: Set transparency

# Overlay data points with bold labels
ax_violin2.scatter(0 + jitter_A, data_A, color=color_A, label="1", alpha=0.7, s=12)
ax_violin2.scatter(0 + jitter_B, data_B, color=color_B, label="2", alpha=0.7, s=12)
ax_violin2.scatter(0 + jitter_C, data_C, color=color_C, label="3", alpha=0.7, s=12)

# Separate points for each violin
ax_violin2.scatter(1 + jitter_A[:4], data_A_1, color=color_A, edgecolor='black', alpha=0.7, marker='s', s=12)
ax_violin2.scatter(1 + jitter_A[4:8], data_A_2, color=color_A, edgecolor='black', alpha=0.7, marker='v', s=12)
ax_violin2.scatter(1 + jitter_A[8:], data_A_3, color=color_A, edgecolor='black', alpha=0.7, marker='*', s=12)

ax_violin2.scatter(2 + jitter_B[:4], data_B_1, color=color_B, edgecolor='black', alpha=0.7, marker='s', s=12)
ax_violin2.scatter(2 + jitter_B[4:8], data_B_2, color=color_B, edgecolor='black', alpha=0.7, marker='v', s=12)
ax_violin2.scatter(2 + jitter_B[8:], data_B_3, color=color_B, edgecolor='black', alpha=0.7, marker='*', s=12)

ax_violin2.scatter(3 + jitter_C[:4], data_C_a_1, color=color_C, edgecolor='black', alpha=0.7, marker='s', label=r'$1^{st}$ cage', s=12)
ax_violin2.scatter(3 + jitter_C[4:8], data_C_a_2, color=color_C, edgecolor='black', alpha=0.7, marker='v', s=12)
ax_violin2.scatter(3 + jitter_C[8:12], data_C_a_3, color=color_C, edgecolor='black', alpha=0.7, marker='*', s=12)
ax_violin2.scatter(3 + jitter_C[12:16], data_C_b_1, color='black', edgecolor=color_C, alpha=0.7, marker='s', label=r'$2^{nd}$ cage', s=12)
ax_violin2.scatter(3 + jitter_C[16:20], data_C_b_2, color='black', edgecolor=color_C, alpha=0.7, marker='v', s=12)
ax_violin2.scatter(3 + jitter_C[20:24], data_C_b_3, color='black', edgecolor=color_C, alpha=0.7, marker='*', s=12)
ax_violin2.scatter(3 + jitter_C[24:28], data_C_c_1, color='white', edgecolor=color_C, alpha=0.7, marker='s', label=r'$3^{rd}$ cage', s=12)
ax_violin2.scatter(3 + jitter_C[28:32], data_C_c_2, color='white', edgecolor=color_C, alpha=0.7, marker='v', s=12)
ax_violin2.scatter(3 + jitter_C[32:], data_C_c_3, color='white', edgecolor=color_C, alpha=0.7, marker='*', s=12)

# Create legend with regular font as default
fig_violin2.canvas.draw()
legend = ax_violin2.legend(
    fontsize=6,
    loc='lower center',
    bbox_to_anchor=(0.5, -0.2),
    ncol=3,
    frameon=False,
    prop=bold_font
)

# Redraw to ensure text objects exist
fig_violin2.canvas.draw()

# Now manually set font properties for each legend text  # update: this doesnt work as the redraw line in the end makes everything regular again. 
# BUT without that line the legend somehow is not boxed. how annoying
legend_texts = legend.get_texts()
for i in range(3):  # First three entries: "1", "2", "3"
    legend_texts[i].set_fontproperties(bold_font)
for i in range(3, len(legend_texts)):  # Remaining entries
    legend_texts[i].set_fontproperties(regular_font)

# Set axis labels
ax_violin2.set_xlabel('Structure')
ax_violin2.set_ylabel('N···N distance (Å)')

# X-tick labels: 'All' (regular), '1/2/3' (bold)
ax_violin2.set_xticks([0,1,2,3])
ax_violin2.set_xticklabels(["All", "1", "2", "3"], fontsize=6)  # normal strings first

# Force a draw so the text objects exist
fig_violin2.canvas.draw()

# Now assign each tick's FontProperties explicitly                # update: not working either, im giving up, making everything regular. this is not going to make it in the paper anyway. 
# for the half violin in the paper I'm using a stupid way, simply printing everything with manual coordinates.
# xtick_labels = ax_violin2.get_xticklabels()
# xtick_labels[0].set_fontproperties(regular_font) #"All" => regular
# xtick_labels[0].set_fontsize(6)
# for lbl in xtick_labels[1:]:
 #    lbl.set_fontproperties(bold_font)             # "1","2","3" => bold
 #    lbl.set_fontsize(6)

# Add legend for the first violin        #update: this line is annoying but necessary
ax_violin2.legend()

# Save the plot
fig_violin2.savefig('violin_plot_with_datapoints.svg', dpi=600, bbox_inches='tight', transparent=True)

# --- Line Plot ---
# Create the second figure for the line plot
fig_line, ax_line = plt.subplots(figsize=(6, 4.5))

# Create the line plot with regular labels
index = np.arange(len(distancesa))
ax_line.plot(index, distancesa, label=r'$\mathit{a}$', color='blue')
ax_line.plot(index, distancesb, label=r'$\mathit{b}$', color='orange')
ax_line.plot(index, distancesc, label=r'$\mathit{c}$', color='green')
ax_line.plot(index, distancesd, label=r'$\mathit{d}$', color='red')

# Add a legend with italic font
ax_line.legend(prop=italic_font)

# Set axis labels for the line plot
ax_line.set_xlabel('Index')
ax_line.set_ylabel('N···N Distance (Å)')

# Add a legend to the line plot
ax_line.legend()

# Save the line plot as a svg file
fig_line.savefig('line_plot_5.svg', dpi=600, bbox_inches='tight')

# plt.close(fig_violin)
plt.close(fig_violin2)
plt.close(fig_line)
