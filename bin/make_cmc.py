#!/usr/bin/env python3
import parmed as pmd

# Create a CMC residue
cmc_residue = pmd.Residue(name='CMC', number=1, chain='A')

# Add atoms to the residue
atom_names = ['C1', 'O1', 'C2', 'O2', 'C3', 'C4', 'O4', 'C5', 'O5', 'C6', 'O6', 'O7', 'HO1', 'HO2', 'HO4', 'HO5', 'HO6']
atom_types = ['C', 'O', 'C', 'O', 'C', 'C', 'O', 'C', 'O', 'C', 'O', 'O', 'H', 'H', 'H', 'H', 'H']

for i in range(len(atom_names)):
    atom = pmd.Atom(name=atom_names[i], type=atom_types[i])
    cmc_residue.add_atom(atom)

# Create a topology object
topology = pmd.Topology()

# Add the CMC residue to the topology
chain = pmd.Chain(name='A')
residue = cmc_residue
topology.add_chain(chain)
chain.add_residue(residue)

# Create a unitcell with 5x5x5 nm dimensions
unitcell = pmd.box_vectors.BoxVectors(5, 5, 5)

# Create a structure object with the topology and unitcell
structure = pmd.Structure()
structure.topology = topology
structure.add_unitcell(unitcell)

# Save the structure to a PDB file
structure.save('cmc.pdb', overwrite=True)
