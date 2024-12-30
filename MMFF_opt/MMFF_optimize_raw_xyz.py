from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem

import os
import glob

# from tqdm import tqdm

def mol_to_xyz(mol, filename):
    with open(filename, 'w') as f:
        num_atoms = mol.GetNumAtoms()
        f.write(f"{num_atoms}\n")
        f.write("MMFF Optimized molecule\n")
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n")

raw_path = "./T1_xyzs"
MMFF_path = "./MMFF_PC_xyz"
MAXITER = 2000

xyz_files = glob.glob(os.path.join(raw_path, '*.xyz'))

for xyz_file in xyz_files:
    if 'W2' in xyz_file:
        os.remove(xyz_file)
        print(xyz_file, ": NO2 - passed")
        continue
    if ('Br' in xyz_file) or ('W10' in xyz_file):
        os.remove(xyz_file)
        print(xyz_file, ": Br - passed")
        continue
    if 'D10' in xyz_file:
        os.remove(xyz_file)
        print(xyz_file, ": D10 - passed")
        continue
    print(xyz_file)
    raw_mol = Chem.MolFromXYZFile(xyz_file)
    mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineBonds(mol,charge=0)

    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        pos = raw_mol.GetConformer().GetAtomPosition(i)
        conformer.SetAtomPosition(i, pos)
    mol.AddConformer(conformer)

    status = AllChem.MMFFOptimizeMolecule(mol, maxIters=MAXITER)
    if status:
        print(status, xyz_file[len(raw_path)+1:])

    save_path = os.path.join(MMFF_path, xyz_file[len(raw_path)+1:])
    mol_to_xyz(mol, save_path)
    
    os.remove(xyz_file)
