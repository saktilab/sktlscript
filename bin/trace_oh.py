#!/usr/bin/env python3
import sys
import numpy as np

with open(sys.argv[1], 'r') as f:
    step = 0
    for line in f:
        # coords = []
        x = []
        y = []
        z = []
        symbols = []
        nat = int(line)
        info = next(f)
        arr = info.split()
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            symbols.append(arr[0])
            x.append(np.array(arr[1], dtype=np.double))
            y.append(np.array(arr[2], dtype=np.double))
            z.append(np.array(arr[3], dtype=np.double))
        active_sites = [int(sys.argv[2])]
        CN = 0
        cutoff = 1.6
        for i in active_sites:
            CN = 0
            for j in range(0, nat):
                if (i == j):
                    continue
                dist = ((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)**0.5
                if (dist <= cutoff):
                    CN += 1
                    # CN = CN + (2*(1-(dist/cutoff)**6)/(1-(dist/cutoff)**12))
            print(CN)