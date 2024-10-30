#!/usr/bin/env python3

from ase.io import read, write
import sys
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.calculators.emt import EMT

atoms = read(sys.argv[1])

atoms.set_calculator(EMT())

MaxwellBoltzmannDistribution(atoms, 300* units.kB)

dyn = VelocityVerlet(atoms, 5 * units.fs)

def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

# Now run the dynamics
printenergy(atoms)
for i in range(20):
    dyn.run(10)
    printenergy(atoms)