#!/usr/bin/env python
import argparse
import os

def create_packmol_input(output_file, box_size, molecules):
    with open(output_file, 'w') as f:
        f.write("tolerance 2.0\n")
        f.write(f"filetype pdb\n")
        f.write(f"output mixture.pdb\n\n")

        for molecule, count in molecules:
            f.write(f"structure {molecule}\n")
            f.write(f"  number {count}\n")
            f.write(f"  inside box 0. 0. 0. {box_size} {box_size} {box_size}\n")
            f.write("end structure\n\n")

def main():
    parser = argparse.ArgumentParser(description="Create a Packmol input file for a mixture of solute molecules in a box.")
    parser.add_argument("-o", "--output", default="packmol_input.inp", help="Output Packmol input file name")
    parser.add_argument("-b", "--box_size", type=float, default=40.0, help="Size of the cubic box")
    parser.add_argument("-m", "--molecules", nargs=2, action='append', metavar=('PDB_FILE', 'COUNT'),
                        help="PDB file and number of molecules (can be used multiple times for different solutes)")

    args = parser.parse_args()

    if not args.molecules:
        parser.error("At least one molecule must be specified.")

    molecules = []
    for pdb_file, count in args.molecules:
        if not os.path.exists(pdb_file):
            parser.error(f"PDB file not found: {pdb_file}")
        try:
            count = int(count)
        except ValueError:
            parser.error(f"Invalid molecule count: {count}")
        molecules.append((pdb_file, count))

    create_packmol_input(args.output, args.box_size, molecules)
    print(f"Packmol input file '{args.output}' has been created successfully.")

if __name__ == "__main__":
    main()
