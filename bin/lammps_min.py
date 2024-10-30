#!/usr/bin/env python3

import sys

from ase.calculators.lammpsrun import LAMMPS


from ase.io import read, write

Atoms = read(sys.argv[1])

calc = LAMMPS()
Atoms.set_calculator(calc)

