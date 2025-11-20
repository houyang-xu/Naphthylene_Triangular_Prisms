'''
This script is used to create conformers of the quadrilateral linker with rotationally flexible arms.

P.C.P. Teeuwen 
Last updated: 11-12-2024
Adjusted for Github: 19-11-2025
'''

from rdkit import Chem
from rdkit.Chem import AllChem
import itertools

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

    # Generate the base molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol, addCoords=True)  # Adds protons to molecule
    AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Creates 3D molecule

    # Extract ranges for independent dihedrals (DA1/DA4 and DA2/DA3 are locked to create a C2 symmetric molecule)
    da1_da4_angles = dihedral_dict["DA1"][4]  # Shared range for DA1 and DA4
    da2_da3_angles = dihedral_dict["DA2"][4]  # Shared range for DA2 and DA3
    angle_combinations = itertools.product(da1_da4_angles, da2_da3_angles) # Create all possible pairs of angles

    # Open the output file for writing all conformers
    with open(output_filename, 'w') as xyz_file:
        for conf_id, (angle1, angle2) in enumerate(angle_combinations):
            conf = mol.GetConformer()
            # Set DA1 and DA4 to angle1
            for key in ["DA1", "DA4"]:
                atom_indices = [idx - 1 for idx in dihedral_dict[key][:4]]  # Adjust to 0-based
                AllChem.SetDihedralDeg(conf, *atom_indices, float(angle1))
            
            # Set DA2 and DA3 to angle2
            for key in ["DA2", "DA3"]:
                atom_indices = [idx - 1 for idx in dihedral_dict[key][:4]]  # Adjust to 0-based
                AllChem.SetDihedralDeg(conf, *atom_indices, float(angle2))

            # Write the resulting conformation to the XYZ file
            xyz_block = Chem.MolToXYZBlock(mol)
            atoms = xyz_block.splitlines()[2:]
            
            xyz_file.write(f"{len(atoms)}\n")  # Number of atoms
            xyz_file.write(f" Conformer {conf_id + 1}: DA1={angle1}, DA2={angle2}\n")  # Comment line
            xyz_file.write("\n".join(atoms) + "\n") 

            print(f"Saved conformer {conf_id + 1} to {output_filename}")

# Example usage
angles_dict = {
    "DA1": [38, 39, 40, 50, range(-90, 270, 10)],  # Atom indices for dihedral 1
    "DA2": [13, 12, 11, 10, range(-90, 270, 10)],  # Atom indices for dihedral 2
    "DA3": [38, 26, 27, 28, range(-90, 270, 10)],  # Atom indices for dihedral 3
    "DA4": [13, 14, 15, 25, range(-90, 270, 10)],  # Atom indices for dihedral 4
}

smiles = "NC1=CC(C=C2)=C(C=C1)C=C2C3=CC(C4=CC=C(C=C(N)C=C5)C5=C4)=C(C6=CC(C=CC(N)=C7)=C7C=C6)C=C3C8=CC=C(C=C(N)C=C9)C9=C8"
smiles_to_xyz_with_locked_dihedrals(smiles, "Conformers_2DA.xyz", angles_dict)

