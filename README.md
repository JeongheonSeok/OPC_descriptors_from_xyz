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

You can also try conformer search during this process (e.g. CREST, https://github.com/grimme-lab/crest), but in my case, it took too much time and I couldn't use it.
If you have enough time, you can try dividing the structure for starting the calculation into neutral/charged particles.

## 1. MMFF optimize
Code for generating suitable initial structures before quantum chemistry calculations. Optimize the structure of .xyz files using MMFF.

1. Create a `raw_xyz` folder and a `MMFF_xyz` folder under the location where `MMFF_optimize_raw_xyz.py` is located.
2. Place .xyz files containing molecular structures to be optimized with MMFF in the `raw_xyz` folder. (Each xyz file contains only one molecule.)
3. run `MMFF_optimize_raw_xyz.py`
```tcsh
python ./MMFF_optimize_raw_xyz.py
```
4. .xyz files containing MMFF optimized structures are stored in `MMFF_xyz` folders.

## 2. Q-Chem Scripts
Code that automates the process of calculating the quantum chemical properties of an OPC from the .xyz files and saving them to csv file.
### 2-1. Calculation Description
The following properties are calculated and saved in the descriptors.csv file generated through this calculation:
1. r_PC: Molecular radius of OPC molecule. Optained from vdW surface area of molecule. $(r\sim \sqrt{\frac{area}{4\pi }})$. in Å
2. Mw: Molecular weight of OPC molecule. in g/mol.
3. E_red: reduction energy. Energy difference between neutral molecule and anion. in eV
4. T1: T1 state energy of neutral OPC molecule. in eV
5. S1: S1 state energy of neutral OPC molecule. in eV
6. S1-T1: Energy difference between S1 state and T1 state of OPC molecule. in eV
7. lambda_int: Internal reorganization energy(of OPC)
8. OS: Oscillator Strength (from TDDFT result)
9. D_T1: The distance between the hole-electron center of the T1 state. $(|\bra r_e - r_h\ket|)$. in Å
10. D_S1: The distance between the hole-electron center of the S1 state. in Å
11. deltaD_T1: 
12. deltaD_S1: 
13. deltaD_(T1-S1): 

The above results are obtained by properly processing Q-Chem out files. If you want to calculate other properties or modify the result output, you can modify `read_rst_in.py`.

The Q-Chem calculation process to obtain these is as follows.

### 2-2. Usage
First, you need to place the scripts and .xyz files for calculation in the appropriate location.
Also before the usage, some of variables in the code must be appropriately modified to suit the user's environment.

1. Place the QCScripts files in the appropriate directory. (This directory must be accessible for future calculations)
2. Copy the file `mk_descriptors_from_xyz.sh` (or `mk_descriptors_from_xyz.pbs` if using the Portable Batch System) to the directory where you will perform the calculation.
3. Create three folders in that directory. Name them `xyz`, `opt`, and `rst`.
4. Put the xyz files to be calculated into the `xyz` folder. (These xyz files have to contain MMFF optimized structure)
5. Modify `mk_descriptors_from_xyz.sh` (or .pbs) appropriately.

   ⑴ In lines 7 and 8, set the QCSCRATCH directory, an environment variable of Q-Chem, to suit the user environment.

   ⑵ In line 16, set the folder directory containing the QCScripts in 1.

   ⑶ If you use PBS, make the same modifications as above for `mk_descriptors_from_xyz.pbs`, but also modify the #PBS part at the top to suit your environment.
6. Run `mk_descriptors_from_xyz.sh`(or .pbs) to perform the calculation.
```tcsh
tcsh ./mk_descriptors_from_xyz.sh
```
or
```tcsh
tcsh ./mk_descriptors_from_xyz.sh > run_output
```
or for pbs usage,
```tcsh
qsub mk_descriptors_from_xyz.pbs
```
7. After the calculation is complete, the results are saved in descriptor.csv.

