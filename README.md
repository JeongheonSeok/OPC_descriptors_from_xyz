# OPC_descriptors_from_xyz
Code for creating a quantum chemistry computation-based descriptor for organic photocatalysts(OPC) (T1, anion-based)

All quantum chemistry computations are done with Q-Chem: 
"Software for the frontiers of quantum chemistry:
   An overview of developments in the Q-Chem 5 package"
   J. Chem. Phys. 155, 084801 (2021)
   https://doi.org/10.1063/5.0055522

## 0. prepare raw xyz structure
In the case of organic photocatalyst molecules, the molecules are often large and complex, so automatic generation of molecular structures using rdkit or open babel was often difficult.
Therefore, the code to automate this process is not uploaded.

You can also try conformer search during this process (e.g. CREST, https://github.com/grimme-lab/crest), but in my case it took too much time and I couldn't use it.
If you have enough time, you can try dividing the structure for starting the calculation into neutral/charged particles.

## 1. Prepare - MMFF optimize
Code for generating suitable initial structures before quantum chemistry calculations. Optimize the structure of .xyz files using MMFF.

1. Create a raw_xyz folder and a MMFF_xyz folder under the location where MMFF_optimize_raw_xyz.py is located.
2. Place .xyz files containing molecular structures to be optimized with MMFF in the raw_xyz folder. (Each xyz file contains only one molecule.)
3. run MMFF_optimize_raw_xyz.py
```bash
python ./MMFF_optimize_raw_xyz.py
```
4. .xyz files containing MMFF optimized structures are stored in MMFF_xyz folders.

## 2. Q-Chem computation
