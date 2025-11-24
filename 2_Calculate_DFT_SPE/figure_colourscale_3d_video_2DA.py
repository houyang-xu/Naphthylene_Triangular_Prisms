'''
This script is used to generate a mp4 video or gif from the conformer data 

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
import matplotlib.animation as animation
import os
import subprocess

# Set font preferences
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica Neue', 'sans-serif']
rcParams['font.size'] = 9  # 9 pt font size

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

# Keep original order - no sorting
dist_a_sequence = dist_a_filtered.values
dist_b_sequence = dist_b_filtered.values
energies_sequence = energies_filtered.values

print(f"Energy range: {energies_sequence.min():.3f} to {energies_sequence.max():.3f} kcal/mol")
print(f"Will animate conformers in original file sequence order")

# Create custom colormap
dark_green = (32/255, 127/255, 76/255)
dark_red = (194/255, 31/255, 48/255)
colors = [dark_green, 'green', 'yellow', 'orange', dark_red]
n_bins = 100
cmap = mcolors.LinearSegmentedColormap.from_list('CustomGreenToRed', colors, N=n_bins)
norm = mcolors.Normalize(vmin=0, vmax=100)

# Create figure with specific size
fig_width = 15 / 2.54
fig_height = 11 / 2.54
#fig_width = 20 / 2.54
#fig_height = 15 / 2.54
fig = plt.figure(figsize=(fig_width, fig_height))
ax = fig.add_subplot(111, projection='3d')

# For large datasets, skip some frames to keep reasonable file size
total_frames = len(energies_sequence)
skip_factor = max(1, total_frames // 900)  # Limit to 900 frames max. actually 811
frame_indices = np.arange(0, total_frames, skip_factor)
print(f"Animation will have {len(frame_indices)} frames (skip factor: {skip_factor})")

# Animation function
def animate(frame_num):
    ax.clear()
    
    # Reset plot properties (they get cleared with ax.clear())
    ax.set_xlim(14.0, 18.0)
    ax.set_ylim(6.25, 12.25)
    ax.set_zlim(0, 100)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.zaxis.set_major_locator(MultipleLocator(20))
    ax.set_xlabel('N···N distance $\it{l}$ (long edge, Å)', fontsize=9, labelpad=5)
    ax.set_ylabel('N···N distance $\it{s}$ (short edge, Å)', fontsize=9, labelpad=5)
    ax.set_zlabel('Relative Energy (kcal·mol⁻¹)', fontsize=9, labelpad=5)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.view_init(elev=20, azim=45)
    
    # Get current frame index (original sequence order)
    frame_idx = frame_indices[frame_num]
    
    # Plot current conformer
    current_x = dist_a_sequence[frame_idx]
    current_y = dist_b_sequence[frame_idx]
    current_z = energies_sequence[frame_idx]
    
    ax.scatter([current_x], [current_y], [current_z], 
              c=[current_z], cmap=cmap, norm=norm, s=50, alpha=0.9, 
              linewidth=1, edgecolors='black')
    
    # Add text annotations
    ax.text2D(0.02, 0.98, f'Energy: {current_z:.2f} kcal/mol', 
              transform=ax.transAxes, fontsize=10, verticalalignment='top',
              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax.text2D(0.02, 0.90, f'Conformer: {frame_idx+1}/{len(energies_sequence)}', 
              transform=ax.transAxes, fontsize=10, verticalalignment='top',
              bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    return ax,

# Create animation
print("Creating animation...")
anim = animation.FuncAnimation(fig, animate, frames=len(frame_indices), 
                              interval=42, blit=False, repeat=True)  # 42ms interval = ~24fps

plt.tight_layout()

# Save as gif
print("Saving as gif file...")
gif_output = 'conformer_animation.gif'
try:
    anim.save(gif_output, writer='pillow', fps=24, dpi=150)
    print(f"gif saved as '{gif_output}'")
    gif_success = True
except Exception as e:
    print(f"gif save fialed: {e}")
    gif_success = False

# Save as mp4
print("Saving as mp4...")
mp4_output = 'conformer_animation.mp4'
mp4_success = False

# Method 1: Try direct ffmpeg
try:
    print("  Trying direct ffmpeg writer...")
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=24, metadata=dict(artist='Conformer Analysis'), bitrate=2400)
    anim.save(mp4_output, writer=writer, dpi=300)
    print(f"mp4 saved as '{mp4_output}' using direct ffmpeg")
    mp4_success = True
except Exception as e:
    print(f"  Direct ffmpeg failed: {e}")

# Method 2: Convert gif to mp4 if gif was successful
if not mp4_success and gif_success:
    try:
        print("  Trying gif to mp4 conversion...")
        cmd = [
            'ffmpeg', '-y', '-i', gif_output, 
            '-movflags', '+faststart', 
            '-pix_fmt', 'yuv420p', 
            '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
            '-r', '24',
            mp4_output
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"mp4 created by converting gif: '{mp4_output}'")
        mp4_success = True
    except subprocess.CalledProcessError as e:
        print(f"  gif to mp4 conversion failed: {e}")
    except FileNotFoundError:
        print("  ffmpeg not found in PATH for conversion")

# Method 3: Save individual frames if needed
if not gif_success and not mp4_success:
    print("both gif and mp4 failed. saving frames only...")
    frames_dir = 'animation_frames'
    os.makedirs(frames_dir, exist_ok=True)
    
    for i, frame_idx in enumerate(frame_indices):
        fig_temp = plt.figure(figsize=(fig_width, fig_height))
        ax_temp = fig_temp.add_subplot(111, projection='3d')
        
        # Set up plot
        ax_temp.set_xlim(14.0, 18.0)
        ax_temp.set_ylim(6.25, 12.25)
        ax_temp.set_zlim(0, 100)
        ax_temp.xaxis.set_major_locator(MultipleLocator(0.5))
        ax_temp.yaxis.set_major_locator(MultipleLocator(1.0))
        ax_temp.zaxis.set_major_locator(MultipleLocator(20))
        ax_temp.set_xlabel('N-N distance $\it{l}$ (long edge, Å)', fontsize=9, labelpad=5)
        ax_temp.set_ylabel('N-N distance $\it{s}$ (short edge, Å)', fontsize=9, labelpad=5)
        ax_temp.set_zlabel('Relative Energy (kcal/mol)', fontsize=9, labelpad=5)
        ax_temp.tick_params(axis='both', which='major', labelsize=8)
        ax_temp.view_init(elev=20, azim=45)
        
        # Plot current conformer
        current_x = dist_a_sequence[frame_idx]
        current_y = dist_b_sequence[frame_idx]
        current_z = energies_sequence[frame_idx]
        
        ax_temp.scatter([current_x], [current_y], [current_z], 
                       c=[current_z], cmap=cmap, norm=norm, s=50, alpha=0.9, 
                       linewidth=1, edgecolors='black')
        
        # Add text
        ax_temp.text2D(0.02, 0.98, f'Energy: {current_z:.2f} kcal/mol', 
                      transform=ax_temp.transAxes, fontsize=10, verticalalignment='top',
                      bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax_temp.text2D(0.02, 0.90, f'Conformer: {frame_idx+1}/{len(energies_sequence)}', 
                      transform=ax_temp.transAxes, fontsize=10, verticalalignment='top',
                      bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(f'{frames_dir}/frame_{i:04d}.png', dpi=300, bbox_inches='tight')
        plt.close(fig_temp)
        
        if (i + 1) % 50 == 0:
            print(f"  Saved {i+1}/{len(frame_indices)} frames")
    
    print(f"All frames saved in '{frames_dir}' directory")
    print("Create video manually with:")
    print(f"  ffmpeg -r 24 -i {frames_dir}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p {mp4_output}")

# Summary
print(f"\nSummary:")
print(f"Total conformers processed: {len(energies_sequence)}")
print(f"Animation frames: {len(frame_indices)}")
print(f"Frame rate: 24 fps")
print(f"Animation duration: {len(frame_indices)/24:.1f} seconds")
print(f"Energy range: {energies_sequence.min():.3f} - {energies_sequence.max():.3f} kcal/mol")
print(f"Sequence: Original file order (no energy sorting)")

if gif_success:
    print(f"gif file: {gif_output}")
if mp4_success:
    print(f"mp4 file: {mp4_output}")

# Don't show the plot window to avoid blocking
plt.close('all')
print("\nAnimation creation complete.")