#!/usr/bin/env python
import subprocess
import sys
import os

def convert_xyz_to_lammps(input_xyz, output_lammps):
    # Convert XYZ to LAMMPS using Open Babel
    subprocess.run(['obabel','-ixyz', input_xyz, '-o', 'lmpdat', '-O', 'temp.lammps', '--nomolecules'])

    # Process the temporary file to remove unwanted sections
    with open('temp.lammps', 'r') as infile, open(output_lammps, 'w') as outfile:
        write_line = True
        for line in infile:
            if any(section in line for section in ['Bonds', 'Angles', 'Dihedrals', 'Impropers']):
                write_line = False
            elif line.strip() == '':
                write_line = True
            if write_line:
                outfile.write(line)

    # Remove the temporary file
    os.remove('temp.lammps')

# Usage
input_file = sys.argv[1]
fname = sys.argv[1].split(".")[0]
output_file = '{}.data'.format(fname)
convert_xyz_to_lammps(input_file, output_file)
