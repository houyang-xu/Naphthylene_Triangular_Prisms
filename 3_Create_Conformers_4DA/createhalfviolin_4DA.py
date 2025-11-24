import numpy as np
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from scipy.stats import gaussian_kde
from matplotlib import rcParams

# Set font preferences with fallback
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
rcParams['font.size'] = 6

# Create font properties for bold and bold italic
bold_font = fm.FontProperties(family='sans-serif', weight='bold')
bold_italic_font = fm.FontProperties(family='Helvetica', weight='bold', style='oblique')

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

            # Write down four N-N distances. This time making them from an rdkit molecule instead of an ase object like in the other scripts. To
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

# Check if the CSV file exists
if os.path.exists('distances_5.csv'):
    # Load the CSV file
    df = pd.read_csv('distances.csv')

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

    # Combine distance arrays into a single 2D array
    data = np.column_stack([distancesa, distancesb, distancesc, distancesd])

    # Convert to DataFrame and save as CSV
    df = pd.DataFrame(data, columns=['DistancesA', 'DistancesB', 'DistancesC', 'DistancesD'])
    df.to_csv('distances_5.csv', index=False)

alldistances = np.concatenate([distancesa, distancesb, distancesc, distancesd])

data = [distancesa, distancesb, distancesc, distancesd]

import matplotlib.pyplot as plt
import numpy as np

# data_A, data_B, data_C from crystals
data_A_1 = [15.4514,10.9164, 15.1883, 10.2954]
data_A_2 = [15.5376, 10.7482, 15.2594, 10.4751]
data_A_3 = [15.3112, 10.9205, 15.2528, 10.3914]
data_A = data_A_1 + data_A_2 + data_A_3

data_B_1 = [16.0221, 9.6617, 15.8152, 9.4057]
data_B_2 = [14.4423, 12.1994, 14.4077, 10.6144]
data_B_3 = [14.886, 12.0502, 14.4883, 10.255]
data_B = data_B_1 + data_B_2 + data_B_3

data_C_a_1 = [16.1086, 9.5407, 15.9601, 9.2241]
data_C_a_2 = [15.8395, 10.0912, 15.6319, 9.832]
data_C_a_3 = [15.8461, 11.072, 15.3791, 9.4627]
data_C_b_1 = [15.9106, 9.9405, 15.7225, 9.6694]
data_C_b_2 = [16.0486, 9.3714, 15.9772, 9.3609]
data_C_b_3 = [16.1358, 9.8791, 15.9161, 9.2164]
data_C_c_1 = [16.2014, 9.4222, 16.1856, 9.1831]
data_C_c_2 = [16.0851, 10.3633, 15.6642, 9.868]
data_C_c_3 = [16.1063, 10.0081, 15.8883, 9.0993]
data_C_a = data_C_a_1 + data_C_a_2 + data_C_a_3
data_C_b = data_C_b_1 + data_C_b_2 + data_C_b_3
data_C_c = data_C_c_1 + data_C_c_2 + data_C_c_3
data_C = data_C_a + data_C_b + data_C_c

# One horizontal half violin, three scatter
fig_hv, ax_hv = plt.subplots(figsize=(8.1768, 2.7))

# 1) Compute KDE for all distances
kde = gaussian_kde(alldistances)
xmin, xmax = np.min(alldistances), np.max(alldistances)
x_vals = np.linspace(xmin, xmax, 300)
pdf_vals = kde(x_vals)

# Function to darken a color (similar to matplotlib's automatic edge darkening)
def darken_color(color, factor=0.7):
    """Darken a color by multiplying RGB values by factor"""
    return tuple(c * factor for c in color)

# Our colours look better
color_A = (0.725, 0.082, 0.227)  # Red tone
color_B = (0.0667, 0.4667, 0.6902)  # Blue tone
color_C = (0, 0.357, 0.188)  # Green tone
color_X = (0.4, 0.4, 0.4)    # Grey tone
edge_color_A = darken_color(color_A)
edge_color_B = darken_color(color_B)
edge_color_C = darken_color(color_C)

# 2) Plot the half violin on top:
#    We fill from y=0 (the "center line") up to y=pdf_vals.
#    The distribution axis is horizontal, so the data is along x.
# ax_hv.fill_between(x_vals, 0, pdf_vals, alpha=0, color=color_X)
# here just to make it narrower
ax_hv.fill_between(x_vals, 0, pdf_vals / 1.5, alpha=1, color=color_X)
# A vertical line at pdf=0 might also be nice, but here we do horizontal approach:
ax_hv.axhline(y=0, color='black', linewidth=1)

# 3) Plot the 3 scatter sets at negative y offsets (below the 0 line).
#    Each scatter is offset more than the previous:
offsetA = -0.05
offsetB = -0.1
offsetC = -0.15
offsetD = 0.035
offsetE = 0.15

# Data go along the x-axis, the offset is along y.

ax_hv.scatter(data_A, np.full_like(data_A, offsetA), color=color_A, edgecolor=edge_color_A, marker='s', s=15, linewidth=1.01, alpha=0.8)
ax_hv.scatter(data_B, np.full_like(data_B, offsetB), color=color_B, edgecolor=edge_color_B, marker='D', s=15, linewidth=1.01, alpha=0.8)
ax_hv.scatter(data_C, np.full_like(data_C, offsetC), color=color_C, edgecolor=edge_color_C, marker = 'o', s=15, linewidth=1.01, alpha=0.8)

# Add bold numbers as labels

ax_hv.text(7, offsetA, "Edges in 1", fontproperties=bold_font, ha="left", va="center", size=6)
ax_hv.text(7, offsetB, "Edges in 2", fontproperties=bold_font, ha="left", va="center", size=6)
ax_hv.text(7, offsetC, "Edges in 3", fontproperties=bold_font, ha="left", va="center", size=6)

# ax_hv.text(8.5, offsetD, "Short Edge (A, C)", fontproperties=bold_font, color="white", ha="left", va="center", size=6)
# ax_hv.text(15.25, offsetD, "Long Edge (B, D)", fontproperties=bold_font, color="white", ha="left", va="center", size=6)

ax_hv.text(8.5, offsetD, "Short Edge (", fontproperties=bold_font, ha="left", va="center", size=6, color="white")
ax_hv.text(9.27, offsetD, "b", fontproperties=bold_italic_font, ha="left", va="center", size=6, color="white")
ax_hv.text(9.36, offsetD, ", ", fontproperties=bold_font, ha="left", va="center", size=6, color="white")
ax_hv.text(9.46, offsetD, "d", fontproperties=bold_italic_font, ha="left", va="center", size=6, color="white")
ax_hv.text(9.56, offsetD, ")", fontproperties=bold_font, ha="left", va="center", size=6, color="white")

ax_hv.text(15.25, offsetD, "Long Edge (", fontproperties=bold_font, ha="left", va="center", size=6, color="white")
ax_hv.text(15.99, offsetD, "a", fontproperties=bold_italic_font, ha="left", va="center", size=6, color="white")
ax_hv.text(16.08, offsetD, ", ", fontproperties=bold_font, ha="left", va="center", size=6, color="white")
ax_hv.text(16.16, offsetD, "c", fontproperties=bold_italic_font, ha="left", va="center", size=6, color="white")
ax_hv.text(16.24, offsetD, ")", fontproperties=bold_font, ha="left", va="center", size=6, color="white")

ax_hv.text(6.5, offsetE, "Occurrence in enumerated conformers", fontproperties=bold_font, color="black", ha="left", va="center", size=6)


# 4) No y labels as meaningless
ax_hv.set_xlim(5.75, 18)
ax_hv.set_ylim(-0.2, 0.2)
ax_hv.set_yticks([])
ax_hv.set_yticklabels([])
ax_hv.set_xlabel("N···N Distance (Å)")

plt.tight_layout()
plt.savefig("half_violin_scatter.svg", dpi=400, bbox_inches='tight', transparent=True)

plt.close(fig_hv)