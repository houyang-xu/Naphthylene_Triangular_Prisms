
# Naphthylene-Triangular-Prisms
This repository contains the scripts and xyz-files relevant for the study of a series of triangular prismatic metal-organic cages (MOCs) from the subcomponent self-assembly of a geometrically flexible quadrilateral tetramine and increasingly large triamine building blocks.

Watch this video to see how the flexible quadrilateral geometry of the tetramine subcomponent works:

<div style="text-align: center; margin: 1em 0; position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%;">
  <iframe style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;" 
          src="https://www.youtube.com/embed/Om998bbISRY"
          frameborder="0" 
          allowfullscreen>
  </iframe>
</div>


## Research Overview
This project investigates the adaptive self-assembly of heteroleptic triangular prismatic cages using a conformationally flexible tetramine subcomponent (**A**) with rotary naphthyl arms. By pairing this adaptive tetramine with three triamines of increasing size (**B**, **C**, **D**), we successfully assembled three distinct triangular prismatic Zn<sub>6</sub>L<sub>3</sub>L'<sub>2</sub> cages with different cavity dimensions and guest-binding properties. The conformational flexibility of the tetramine enables it to adjust its quadrilateral geometry to match triamines of varying sizes, overcoming traditional size-matching constraints in heteroleptic self-assembly.

## Software Used
- [Python](https://www.python.org/) 3.12.8
- [RDKit](https://www.rdkit.org/) for molecular structure generation and manipulation
- [ORCA](https://www.faccts.de/orca/) v6.1.0 for DFT calculations
- [NumPy](https://numpy.org/) for numerical computations
- [Pandas](https://pandas.pydata.org/) for data handling
- [Matplotlib](https://matplotlib.org/) for visualisation
- [ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/) for trajectory analysis

For installing these tools, refer to their respective documentation or package managers suitable for your system.

## Repository Structure
- `scripts/`: Contains all Python scripts for conformational analysis and data visualisation
- `conformers/`: Contains CSV files with conformer geometry data
  - `distances_5.csv`: Complete conformer enumeration (72⁴ = 26,873,856 conformers)
  - `distances_unique_*.csv`: Filtered conformer datasets
  - `conformer_distances_data.csv`: Conformer data with DFT energies
- `DFT_calculations/`: Contains input and output files for DFT calculations
  - *C*<sub>2</sub>-symmetric conformers (1,296 conformers with single-point DFT energies)
  - Non-*C*<sub>2</sub> conformers (filtered by 100 kcal mol⁻¹ threshold)

## Cages
The three triangular prismatic cages (**1**-**3**) were synthesised via subcomponent self-assembly:
- Cage **1**: Assembled from tetramine **A** + triamine **B** (smallest)
- Cage **2**: Assembled from tetramine **A** + triamine **C** (intermediate)
- Cage **3**: Assembled from tetramine **A** + triamine **D** (largest)

## Conformational Analysis

### 1. *C*<sub>2</sub>-Symmetric Conformer Enumeration
Based on previous knowledge that heteroleptic triangular prisms feature *C*<sub>2</sub>-symmetric ligands, an initial conformer search was performed with *C*<sub>2</sub> symmetry constraints.

**Method:**
- Dihedral angles DA1 and DA4, DA2 and DA3 constrained as identical pairs
- Each dihedral angle varied from -90° to 260° in 10° increments
- Total conformers: 36<sup>2</sup> = 1,296

**Computational Details:**
```
Input: SMILES string of tetramine A
Method: ETKDG (RDKit) for 3D structure generation
DFT level: PBE-D3/def2-TZVP with TightSCF convergence
Energy threshold: 100 kcal mol<sup>-1</sup> above minimum
Conformers analysed: 811 (after filtering)
```

**Results:**
- Identified lowest energy conformer (#556: DA1 = DA2 = DA3 = DA4 = 60°)
- 811 conformers within 100 kcal mol<sup>-1</sup> of minimum energy
- Results visualised in scatter plots showing energy-geometry correlation

### 2. Non-*C*<sub>2</sub> Conformer Enumeration
Following discovery of non-*C*<sub>2</sub> conformations in cage **2**'s crystal structure, expanded conformational analysis was performed without symmetry constraints.

**Method:**
- All four dihedral angles (DA1-DA4) varied independently
- Each dihedral angle varied from 0° to 355° in 5° increments
- Total conformers: 72<sup>4</sup> = 26,873,856

**Computational Details:**
```
Method: Systematic dihedral angle enumeration
Increment: 5° per dihedral angle
N-N distance calculation: Direct from Cartesian coordinates
Output: distances_5.csv (2 GB file)
```

**Geometric Parameters:**
For each conformer, four N-N distances define the quadrilateral geometry:
- Distance **A**: N<sub>1</sub>-N<sub>2</sub> (long edge)
- Distance **B**: N<sub>2</sub>-N<sub>3</sub> (short edge)
- Distance **C**: N<sub>3</sub>-N<sub>4</sub> (long edge)
- Distance **D**: N<sub>4</sub>-N<sub>1</sub> (short edge)
- Average long edge: *l* = (**A** + **C**)/2
- Average short edge: *s* = (**B** + **D**)/2

**Results:**
- Visualised as violin plots showing N-N distance distributions
- Scatter plots showing *s* vs *l* reveal conformational landscape
- Identified wide range of accessible geometries (9.3-16.2 Å)

### 3. Conformer Analysis with Single-Crystal X-Ray Diffraction (SCXRD) structures
Analysis of tetramine **A** geometries observed in crystal structures of cages **1**-**3**:

**Cage 1:**
- All three tetramine units in portrait orientation
- Mean geometry: (10.62 ± 0.27) × (15.33 ± 0.13) Å
- More equilateral than default conformation

**Cage 2:**
- Mixed orientations: one portrait, two landscape
- Heterochiral metal vertices (2:1 or 1:2 ratio of Λ:Δ)
- Landscape-oriented quadrilateral ligands: (15.92 ± 0.10) × (9.53 ± 0.13) Å
- Portrait-oriented quadrilateral ligands: (14.56 ± 0.22) × (12.12 ± 0.11, 10.43 ± 0.18) Å
- Most distorted prismatic geometry

**Cage 3:**
- All three tetramine units in landscape orientation
- Mean geometry: (15.92 ± 0.22) × (9.70 ± 0.50) Å
- Close to default conformation
