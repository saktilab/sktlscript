#!/usr/bin/env python
import sys

def cartesian_to_lammps_mol(input_file, output_file):
    # Read input file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Parse input data
    num_atoms = int(lines[0].strip())
    atoms = []
    for line in lines[2:]:  # Skip the first two lines
        parts = line.split()
        if len(parts) == 4:
            atom_type = parts[0]
            x, y, z = map(float, parts[1:])
            atoms.append((atom_type, x, y, z))

    # Write LAMMPS molecule file
    with open(output_file, 'w') as f:
        f.write(f"# LAMMPS molecule file converted from {input_file}\n\n")
        f.write(f"{num_atoms} atoms\n")
        f.write("0 bonds\n")
        f.write("0 angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")

        f.write("Coords\n\n")
        for i, (_, x, y, z) in enumerate(atoms, 1):
            f.write(f"{i} {x:.10f} {y:.10f} {z:.10f}\n")

        f.write("\nTypes\n\n")
        atom_types = {atom_type: i for i, atom_type in enumerate(set(atom[0] for atom in atoms), 1)}
        for i, (atom_type, _, _, _) in enumerate(atoms, 1):
            f.write(f"{i} {atom_types[atom_type]}\n")


input_file = sys.argv[1]
fname = input_file.split(".")[0]
output_file = "{}.mol".format(fname)
cartesian_to_lammps_mol(input_file, output_file)
print(f"Conversion complete. LAMMPS molecule file saved as {output_file}")
