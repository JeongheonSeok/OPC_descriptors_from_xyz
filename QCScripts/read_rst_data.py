import csv
import numpy as np

# Set out file path
qchem_out_file_path = './rst/PCnamehere_rst.out'
csvname = "descriptor.csv"

h2ev = 27.211386246

atomic_weights = {
    "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811, "C": 12.0107,
    "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797, "Na": 22.98977, "Mg": 24.305,
    "Al": 26.981538, "Si": 28.0855, "P": 30.973761, "S": 32.065, "Cl": 35.453, "K": 39.0983,
    "Ar": 39.948, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Ni": 58.693, "Co": 58.933, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95,
    "Tc": 98, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "I": 126.90447, "Te": 127.60, "Xe": 131.29,
    "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
    "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.04, "Lu": 174.97, "Hf": 178.49,
    "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Th": 232.04,
    "Pa": 231.04, "U": 238.03, "Np": 237, "Pu": 244, "Am": 243, "Cm": 247, "Bk": 247,
    "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259, "Lr": 262, "Rf": 267,
    "Db": 270, "Sg": 271, "Bh": 270, "Hs": 277, "Mt": 276, "Ds": 281, "Rg": 282,
    "Cn": 285, "Nh": 286, "Fl": 289, "Mc": 290, "Lv": 293, "Ts": 294, "Og": 294
}       # See python periodictable library
PI = 3.14159265358979

def read_qchem_output(file_path):
    # Read 6 Total energy, 'T1', 'S1' from the top of the .out file of file_path and return them as total_energies array.
    total_energies = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        total_energy_count = 0
        for line in lines:
            if 'Total energy' in line:
                total_energy = float(line.split('=')[1])
                total_energies.append(total_energy)
                total_energy_count += 1
                if total_energy_count == 6:
                    break
        total_energies = [total_energies[5]]+total_energies[:5]   # The changed rst.in format omits the neutral state SP calculation and obtains the corresponding information from the TDDFT calculation. The order of the list is changed because the value that should have been read first is read sixth.
        total_energies.append(0);     total_energies.append(0)
        for line in lines:
            if 'Excited state' in line and 'excitation energy (eV)' in line:
                excitation_energy = float(line.split('=')[1])
            if 'Multiplicity: Triplet' in line:
                total_energies[6] = excitation_energy
            if 'Multiplicity: Singlet' in line:
                total_energies[7] = excitation_energy
            if (total_energies[6]!=0) & (total_energies[7]!=0):
                break
        for line in lines:
            if 'Excited state   3' in line:
                print("There are more excited state than 2")
    return total_energies

def calculate_molecular_weight(file_path):
    # Calculate the molecular weight
    xyzfiledir = './xyz/'+file_path[6:-8]+'.xyz'
    with open(xyzfiledir, 'r') as file:
        lines = file.readlines()
    num_atoms = int(lines[0].strip())
    molecular_weight = 0.0
    for i in range(2, 2 + num_atoms):
        atom_type = lines[i].split()[0]
        molecular_weight += atomic_weights[atom_type]
    return molecular_weight

def calc_rPC(file_path):
    with open(file_path, 'r') as file:
        count = False
        for line in file:
            if "Molecular Surface Area" in line:
                if count:
                    molecular_surface_area = float(line.split('=')[1].strip().split()[0])
                    break
                else:
                    count = True
    if molecular_surface_area is not None:
        r_PC = (molecular_surface_area/(4*PI))**0.5
    else: r_PC = 0
    return r_PC

def read_OS(file_path):
    strength = 0.
    with open(file_path, 'r') as file:
        count = False
        for line in file:
            if count & ('Strength' in line):
                strength = float(line.split(':')[1].strip().split()[0])
                break
            if 'Multiplicity: Singlet' in line:
                count = True
    return strength

def read_D(file_path):
    Ds = [-1., -1.]
    with open(file_path, 'r') as file:
        AUD_check = False
        ST_check = 'yet'
        for line in file:
            if "Analysis of Unrelaxed Density Matrices" in line:
                AUD_check = True
            if AUD_check:
                if (Ds[0]!=-1.) & (Ds[1]!=-1.):
                    break
                if 'Triplet' in line:
                    ST_check = 'T'
                if 'Singlet' in line:
                    ST_check = 'S'
                if ST_check == 'T':
                    if '|<r_e - r_h>|' in line:
                        Ds[0] = float(line.split(':')[1].strip().split()[0])
                        ST_check = 'yet'
                if ST_check == 'S':
                    if '|<r_e - r_h>|' in line:
                        Ds[1] = float(line.split(':')[1].strip().split()[0])
                        ST_check = 'yet'
    return Ds

def read_deltaD(file_path):
    deltaDs = [0., 0., 0.]
    with open(file_path, 'r') as file:
        ESA_check = False
        D_GS = [0.,0.,0.];  D_T1 = -1;  D_S1 = -1
        AUD_check = False
        ST_check = 'yet'
        for line in file:
            if 'Excited State Analysis' in line:
                ESA_check = True
            if ESA_check:
                if 'Cartesian components [D]:' in line:
                    temp = line.split(':')[1].strip().split()
                    if len(temp) == 4:
                        temp.pop(0)
                    elif len(temp) == 3:
                        temp[0] = temp[0][1:]
                    D_GS = [float(temp[0][:-1]), float(temp[1][:-1]), float(temp[2][:-1])]
                    ESA_check = False
            if 'Analysis of Unrelaxed Density Matrices' in line:
                AUD_check = True
            if AUD_check:
                if (D_T1!=-1) & (D_S1!=-1):
                    break
                if 'Triplet' in line:
                    ST_check = 'T'
                if 'Singlet' in line:
                    ST_check = 'S'
                if ST_check == 'T':
                    if 'Cartesian components [D]:' in line:
                        temp = line.split(':')[1].strip().split()
                        if len(temp) == 4:
                            temp.pop(0)
                        elif len(temp) == 3:
                            temp[0] = temp[0][1:]
                        D_T1 = [float(temp[0][:-1]), float(temp[1][:-1]), float(temp[2][:-1])]
                        ST_check = 'yet'
                if ST_check == 'S':
                    if 'Cartesian components [D]:' in line:
                        temp = line.split(':')[1].strip().split()
                        if len(temp) == 4:
                            temp.pop(0)
                        elif len(temp) == 3:
                            temp[0] = temp[0][1:]
                        D_S1 = [float(temp[0][:-1]), float(temp[1][:-1]), float(temp[2][:-1])]
                        ST_check = 'yet'
    deltaDs[0] = float(np.linalg.norm(np.array(D_T1)-np.array(D_GS)))
    deltaDs[1] = float(np.linalg.norm(np.array(D_S1)-np.array(D_GS)))
    deltaDs[2] = float(np.linalg.norm(np.array(D_T1)-np.array(D_S1)))
    return deltaDs

def calc_lambda_internal(qcout_E):
    lambda_internal = 0.5*(qcout_E[2]+qcout_E[3]-qcout_E[0]-qcout_E[1])*h2ev
    return lambda_internal

# Read total energy and excitation energy from Q-Chem out file
qcout_E = read_qchem_output(qchem_out_file_path)
energies = [(qcout_E[5]-qcout_E[4])*h2ev, qcout_E[6], qcout_E[7], qcout_E[7]-qcout_E[6]]
energies.append(calc_lambda_internal(qcout_E))
energies.append(read_OS(qchem_out_file_path))
Ds = read_D(qchem_out_file_path)
energies.append(Ds[0]); energies.append(Ds[1])
deltaDs = read_deltaD(qchem_out_file_path)
energies.append(deltaDs[0]);    energies.append(deltaDs[1]);    energies.append(deltaDs[2])
energies.insert(0, calculate_molecular_weight(qchem_out_file_path))
energies.insert(0, calc_rPC(qchem_out_file_path))
energies.insert(0, qchem_out_file_path[6:-8])

with open(csvname, 'a', newline='') as csvfile:
    # fieldnames = ["PC name", "r_PC", "Mw", "E_nn", "E_aa", "E_an", "E_na", "E_n_full", "E_a_full", "T1", "S1"]
    fieldnames = ["PC name", "r_PC", "Mw", "E_red", "T1", "S1", "S1-T1", "lambda_int", "OS", "D_T1", "D_S1", "deltaD_T1", "deltaD_S1", "deltaD_(T1-S1)"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write header if file is empty
    if csvfile.tell() == 0:
        writer.writeheader()

    # Add data
    writer.writerow({fieldnames[i]: energies[i] for i in range(len(energies))})
