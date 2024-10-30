#!/usr/bin/env python3

import sys


spin_polarized = False
SPIN_UNKNOWN = 0
SPIN_UP = 1
SPIN_DOWN = 2

current_spin = SPIN_UNKNOWN

energy_level = []
energy_level_alpha = []
energy_level_beta = []

filling = []
filling_alpha = []
filling_beta = []

infile = 'detailed.out'

scc = True

if (len(sys.argv) > 1):
    infile = sys.argv[1]

with open(infile, 'r') as f:
    for line in f:
        if 'COMPONENT' in line:
            arr = line.strip().split(' = ')
            if arr[1] == 'q':
                spin_polarized = False
            elif arr[1].lower() == 'up':
                spin_polarized = True
                current_spin = SPIN_UP
                energy_level = energy_level_alpha
                filling = filling_alpha
            elif arr[1].lower() == 'down':
                spin_polarized = True
                current_spin = SPIN_DOWN
                energy_level = energy_level_beta
                filling = filling_beta
            else:
                raise ValueError('COMPONENT value is wrong')
        if 'Eigenvalues /H' in line:
            while (True):
                line = next(f).strip()
                if not line:
                    break
                value = float(line)
                energy_level.append(value)
        if 'Fillings' in line:
            while (True):
                line = next(f).strip()
                if not line:
                    break
                value = float(line)
                filling.append(value)
        if 'SCC is NOT converged' in line:
            scc = False

if (spin_polarized):
    for ev_a, fil_a, ev_b, fil_b in zip(energy_level_alpha, filling_alpha, energy_level_beta, filling_beta):
        print ('{:20.12f} {:10.5f} {:20.12f} {:10.5f}'.format(ev_a, fil_a, ev_b, fil_b))
else:
    for ev, fil in zip(energy_level, filling):
        print ('{:20.12f} {:10.5f}'.format(ev, fil))

if (not scc):
    print ('WARNING: SCC NOT CONVERGED')
