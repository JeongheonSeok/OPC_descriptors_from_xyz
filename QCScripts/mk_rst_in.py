import os
import re

def parse_qchem_output(file_content):
    """Extract optimized nuclear coordinates from Q-Chem output."""
    converged_pattern = re.compile(r'\*\*  OPTIMIZATION CONVERGED  \*\*')
    orientation_pattern = re.compile(r'Standard Nuclear Orientation \(Angstroms\)')
    lines = file_content.split('\n')
    optimized_sections = []
    
    i = 0
    while i < len(lines):
        if converged_pattern.search(lines[i]):
            # Find the start of the Standard Nuclear Orientation section
            while i < len(lines) and not orientation_pattern.search(lines[i]):
                i += 1
            if i < len(lines):
                start_idx = i + 3  # Skip header lines
                i = start_idx
                # Find the end of the coordinates section
                section = []
                while i < len(lines) and not lines[i].startswith(' ----------------------------------------------------------------'):
                    if lines[i].strip():  # Ignore empty lines
                        section.append(' '.join(lines[i].split()[1:]))
                    i += 1
                optimized_sections.append(section)
        i += 1
    
    return optimized_sections

def create_in_file(template_content, neu_coord, ani_coord, output_filename):
    """Create a .in file by replacing placeholders with coordinates."""
    neu_coord_str = "\n".join(neu_coord)
    ani_coord_str = "\n".join(ani_coord)
    
    new_content = template_content.replace("neu_xyz", neu_coord_str).replace("ani_xyz", ani_coord_str)
    
    with open(output_filename, 'w') as f:
        f.write(new_content)

opt_file_path = './opt/PCnamehere_opt.out'
rst_path = './rst/'
os.makedirs(rst_path, exist_ok=True)
template_path = './rst_format.in'
output_filename = './rst/PCnamehere_rst.in'

with open(template_path, 'r') as f:
    template_content = f.read()
with open(opt_file_path, 'r') as f:
    file_content = f.read()
coordinate_sections = parse_qchem_output(file_content)
if len(coordinate_sections) >= 2:
    neu_coord = coordinate_sections[0]
    ani_coord = coordinate_sections[1]
    create_in_file(template_content, neu_coord, ani_coord, output_filename)
else:
    print("File PCnamehere_opt.out does not contain enough optimized sections.")
