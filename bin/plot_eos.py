#!/usr/bin/env python3
from ase.eos import *
import sys
import argparse
import os
from ase.units import kJ
from ase.io import write,read 
from ase import Atoms

parser = argparse.ArgumentParser(description='Fit the E-V curve to equation of state (EOS)')
parser.add_argument('-i', '--input', type=str, help='File containing two columns of volume and energy, respectively')
parser.add_argument('-eos', '--eos', type=str, help='Type of EOS, i.e., sj, taylor, murnaghan, birchmurnaghan, pouriertarantola, vinet, antonschmidt, and p3')
opt = parser.parse_args(sys.argv[1:])

f = open(opt.input, 'r')
volumes = []
energies = []
evperAngtoGPa = 160.21766208
hartree2ev = 27.2114
natom = f.readline()
for line in f:
    arr=line.split()
    if (os.path.exists("./dftb_in.hsd")):
        energy = float(arr[1])*hartree2ev/float(natom)
    elif (os.path.exists("./INCAR")):
        energy = float(arr[1])/float(natom)
    volumes.append(float(arr[0]))
    energies.append(energy)

eos = EquationOfState(volumes, energies, eos=opt.eos)
v0, e0, B = eos.fit()

print('B = {} GPa'.format(B / kJ * 1.0e24))
# print('a = {} Angstrom'.format(a))
eos.plot("eos.png")

