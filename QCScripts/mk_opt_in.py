import os

################################User Settings#################################
# There should be Q-Chem input files, PBS files, and xyz files in the folder.

directory = '.'             # Set current directory path
xyzdir = './xyz'
optdir = './opt'
qcformat = 'opt_format.in'
xyzindicator = 'xyzloc'             # Replace the text in the input template with coordinates
##############################################################################

# Load .xyz file name
xyz_files = []
for filename in os.listdir(xyzdir):
    if filename.endswith('.xyz'):
        xyz_files.append(filename)


# Read the contents of QC format and save them to qcf_txt.
with open(qcformat, 'r') as qcf:
    qcf_txt = qcf.read()

for xyz_file in xyz_files:
    # Read coordinate information from xyz file and save it in coord.
    with open(os.path.join(xyzdir, xyz_file), 'r') as xf:
        lines_xf = xf.readlines()
    coord = ''.join(lines_xf[2:])
        
    # Assign coord to format.
    qcin = qcf_txt
    qcin = qcin.replace(xyzindicator, coord)
        
    # Create an input file with coordinates entered
    file_name = xyz_file[:-4] + '_opt.in'
    if not os.path.exists(optdir):
        os.makedirs(optdir)
    folder_path = xyz_file[:-4] + '/opt'
    file_path = os.path.join(optdir, file_name)
    with open(file_path, 'w') as file:
        file.write(qcin)
