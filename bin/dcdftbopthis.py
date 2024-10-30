#!/usr/bin/env python

import sys

logfile = sys.argv[1]

nat = 0
step = 1
with open(logfile, 'r') as f:
    for line in f:
        if 'Total number of atoms' in line:
            arr = line.split()
            nat = int(arr[5])
        if 'Molecular coordinate [Angstrom]' in line:
            next(f)
            next(f)
            next(f)
            next(f)
            print(nat)
            print('Optimization Step:{}'.format(step))
            step += 1
            for i in range(nat):
                line = next(f)
                print(line.strip())

