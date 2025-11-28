'''
This script uses "distances_5.csv" created by "createconformers_4DA.py" to create a scatter plot, with the over 2.6*10^8 conformers enumerated represented with the
mean long edge (l), and mean short edge (s), as a method to show the overall theoretically accessible quadrilateral geometry of tetramine A. In comparison with the 
Nimine···Nimine Distances (Å) in cages, showing the broad scope of quadrilateral geometries that A adopts to show the geometric versatility of A.

P.C.P. Teeuwen and H. Xu
Lasted updated by Houyang Xu
Last updated: 1 Aug 2025
Updated for GitHub: 24 Nov 2025
'''

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

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

# Use the regular font as the default for all text (avoid math text formatting)
mpl.rcParams['mathtext.fontset']    = 'custom'
mpl.rcParams['mathtext.rm']         = regular_font.get_name()
mpl.rcParams['mathtext.bf']         = bold_font.get_name()
mpl.rcParams['mathtext.default']    = 'regular'
mpl.rcParams['text.usetex']         = False
mpl.rcParams['font.family']         = regular_font.get_name()
mpl.rcParams['font.size']           = 6

# Read the precomputed points
df = pd.read_csv('distances_5.csv')
distancesa = df['DistancesA'].values
distancesb = df['DistancesB'].values
distancesc = df['DistancesC'].values
distancesd = df['DistancesD'].values
print("Data loaded successfully.")

# Calculate average lengths l and s here. For my own version I made another .csv file with 2 columns of averaged l and s values so theres no need to calculate
# average every time making the plot.
x_gray = (distancesa + distancesc) / 2.0
y_gray = (distancesb + distancesd) / 2.0

# Another .csv file with pre-calculated average l and s values
# df_gray = pd.read_csv('grey_points.csv')
# x_gray = df_gray['Xgray'].values
# y_gray = df_gray['Ygray'].values

# Scatters here

data_A_1 = [15.4514, 10.9164, 15.1883, 10.2954]
data_A_2 = [15.5376, 10.7482, 15.2594, 10.4751]
data_A_3 = [15.3112, 10.9205, 15.2528, 10.3914]
data_A   = [data_A_1, data_A_2, data_A_3]

data_B_1 = [16.0221,  9.6617, 15.8152,  9.4057]
data_B_2 = [14.4423, 12.1994, 14.4077, 10.6144]
data_B_3 = [14.886,  12.0502, 14.4883, 10.255]
data_B   = [data_B_1, data_B_2, data_B_3]

data_C_a_1 = [16.1086,  9.5407, 15.9601,  9.2241]
data_C_a_2 = [15.8395, 10.0912, 15.6319,  9.832]
data_C_a_3 = [15.8461, 11.072,  15.3791,  9.4627]
data_C_b_1 = [15.9106,  9.9405, 15.7225,  9.6694]
data_C_b_2 = [16.0486,  9.3714, 15.9772,  9.3609]
data_C_b_3 = [16.1358,  9.8791, 15.9161,  9.2164]
data_C_c_1 = [16.2014,  9.4222, 16.1856,  9.1831]
data_C_c_2 = [16.0851, 10.3633, 15.6642,  9.868]
data_C_c_3 = [16.1063, 10.0081, 15.8883,  9.0993]
data_C = [
    data_C_a_1, data_C_a_2, data_C_a_3,
    data_C_b_1, data_C_b_2, data_C_b_3,
    data_C_c_1, data_C_c_2, data_C_c_3
]

def avg_first_third_second_fourth(arr4):
    return ((arr4[0] + arr4[2]) / 2.0, (arr4[1] + arr4[3]) / 2.0)

points_A = [avg_first_third_second_fourth(arr) for arr in data_A]
points_B = [avg_first_third_second_fourth(arr) for arr in data_B]
points_C = [avg_first_third_second_fourth(arr) for arr in data_C]

#plt.figure(figsize=(7,5))
plt.figure(figsize=(3.7402,2.7559))

# fig_width = 17 / 2.54
# fig_height = 12.75 / 2.54
# plt.figure(figsize=(fig_width, fig_height))

# Gray points (label "Theoretical" remains in regular font)
plt.scatter(x_gray, y_gray, color='gray', alpha=0.5, s=1, rasterized=True, label='Theoretical')

# Define colours
color_A = (0.725, 0.082, 0.227)  # Red tone
color_B = (0.0667, 0.4667, 0.6902)  # Blue tone
color_C = (0, 0.357, 0.188)  # Green tone

# Function to darken a color (similar to matplotlib's automatic edge darkening)
def darken_color(color, factor=0.7):
    """Darken a color by multiplying RGB values by factor"""
    return tuple(c * factor for c in color)

# Get darker versions for edges
edge_color_A = darken_color(color_A)
edge_color_B = darken_color(color_B)
edge_color_C = darken_color(color_C)

# Plot the points with legend labels (plain text for numbers)
for i, (xA, yA) in enumerate(points_A):
    lbl = "1" if i == 0 else None
 #    lbl = "4.11" if i == 0 else None
  #  plt.scatter(xA, yA, color=color_A, edgecolor=edge_color_A, marker='s', s=36, label=lbl)
    plt.scatter(xA, yA, color=color_A, edgecolor=edge_color_A, marker='s', s=15, label=lbl)

for i, (xB, yB) in enumerate(points_B):
    lbl = "2" if i == 0 else None
    #    lbl = "4.12" if i == 0 else None
   #   plt.scatter(xB, yB, color=color_B, edgecolor=edge_color_B, s=36, label=lbl)
    plt.scatter(xB, yB, color=color_B, edgecolor=edge_color_B, s=15, label=lbl)

for i, (xC, yC) in enumerate(points_C):
    lbl = "3" if i == 0 else None
   #  lbl = "4.13" if i == 0 else None
   #    plt.scatter(xC, yC, color=color_C, edgecolor=edge_color_C, marker='^', s=36, label=lbl)
    plt.scatter(xC, yC, color=color_C, edgecolor=edge_color_C, marker='^', s=15, label=lbl)

# Use plain text for axis labels (avoiding math text) and set font size and font explicitly
# plt.xlabel("N–N distance a (long edge, Å)", fontsize=6, fontproperties=regular_font)
# plt.ylabel("N–N distance b (long edge, Å)", fontsize=6, fontproperties=regular_font)

#plt.xlabel(r"N–N distance $\mathit{a}$ (long edge, Å)", fontsize=12, fontproperties=regular_font)
#plt.ylabel(r"N–N distance $\mathit{b}$ (short edge, Å)", fontsize=12, fontproperties=regular_font)

# Use plain text for axis labels (avoiding math text) and set font size and font explicitly
plt.xlabel(r"N···N distance $\mathit{l}$ (long edge, Å)", fontsize=6, fontproperties=regular_font)
plt.ylabel(r"N···N distance $\mathit{s}$ (short edge, Å)", fontsize=6, fontproperties=regular_font)

# plt.xlim(13.25, 17.75)
# plt.ylim(5.75, 12.25)

# ticks at every integer between the limits (inclusive)
# plt.xticks(np.arange(14, 18, 1), fontsize=6, fontproperties=regular_font)  # 14, 15, 16, 17
# plt.yticks(np.arange(6, 13, 1), fontsize=6, fontproperties=regular_font)   # 6, 7, …, 12
plt.xticks(fontsize=6, fontproperties=regular_font)
plt.yticks(fontsize=6, fontproperties=regular_font)

# Create legend with default regular font
leg = plt.legend(loc='best', fontsize=6, prop=regular_font)

# Update legend texts so that "1", "2", "3" are bold while "Theoretical" remains regular
for text in leg.get_texts():
    text.set_fontsize(6)
    if text.get_text() in ["1", "2", "3"]:
 #   if text.get_text() in ["4.11", "4.12", "4.13"]:
        text.set_fontproperties(bold_font)
    else:
        text.set_fontproperties(regular_font)

plt.tight_layout()
# plt.savefig('scatter_AC_vs_BD_pres.svg', dpi=400, transparent=True)
plt.savefig('scatter_AC_vs_BD_4DA.svg', dpi=400, transparent=True)
plt.close()