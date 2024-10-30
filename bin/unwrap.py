#!/usr/bin/env python3

import MDAnalysis as mda
import sys

# Load the XYZ trajectory
u = mda.Universe(sys.argv[1])

# Unwrap the trajectory
u.atoms.wrap(compound='segments', inplace=True)

# Save the unwrapped trajectory
u.atoms.write("output_unwrapped.xyz")

