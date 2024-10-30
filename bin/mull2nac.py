#!/usr/bin/env python3

import sys
import numpy as np

file = sys.argv[1]
fout = sys.stdout
with open(file, 'r') as fin:
    for line in fin:
        norb, nat = list(map(int,line.split()))
        dummy = next(fin)
        arr = dummy.split()

        symbols = {}
        nac = np.zeros(nat, np.double)
        for j in range(norb):
            line = next(fin)
            arr = line.split()
            atom_index = int(arr[0])-1
            nac[atom_index] += float(arr[3])
            symbols[atom_index] = arr[1]
        
        print(len(nac), file=fout)
        print(dummy, end='', file=fout)
        for i in range(len(nac)):
            print('{:9d} {:20.12f}'.format(i+1, nac[i]), file=fout)

