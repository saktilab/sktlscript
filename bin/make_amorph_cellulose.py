#!/usr/bin/env python3
import openbabel as ob

# Set up molecule
mol = ob.OBMol()

# Add cellulose monomer to molecule
glucose = ob.OBMol()
obConversion = ob.OBConversion()
obConversion.SetInFormat('mol')
obConversion.ReadFile(glucose, 'glucose.mol')
mol += glucose

# Copy and translate monomer to generate amorphous structure
num_monomers = 1000  # number of monomers to include
displacement = 0.5  # maximum displacement for each monomer
for i in range(num_monomers):
    # Copy monomer
    new_mol = ob.OBMol(mol)

    # Randomly displace each atom
    for atom in ob.OBMolAtomIter(new_mol):
        x, y, z = atom.GetVector()
        x += (2*ob.random() - 1) * displacement
        y += (2*ob.random() - 1) * displacement
        z += (2*ob.random() - 1) * displacement
        atom.SetVector(x, y, z)

    # Add copied and displaced monomer to original molecule
    mol += new_mol

# Save molecule to file
obConversion.SetOutFormat('xyz')
obConversion.WriteFile(mol, 'amorphous_cellulose.xyz')
