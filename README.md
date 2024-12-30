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

If you have the time, you can try dividing the structure for starting the calculation into neutral/charged particles.

## 1. MMFF optimize
