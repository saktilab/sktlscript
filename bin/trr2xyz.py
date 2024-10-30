#!/usr/bin/env python3
import MDAnalysis as mda

pdb = "run.pdb"

u = mda.Universe(pdb)
with mda.Writer("run.xyz", n_atoms=u.atoms.n_atoms) as xyz:
    for ts in u.trajectory:
        xyz.write(ts)
