#!/usr/bin/env python3
import sys
import numpy as np

E_rep = []
r = []
with open(sys.argv[1], "r") as f:
    for line in f:
        arr = line.split()
        E_rep.append(arr[1])
        r.append(arr[0])
    E_rep.append(0.0)    
    r.append(float(sys.argv[2]))
    E_rep = np.array(E_rep, dtype='float64')
    
    r = np.array(r, dtype='float64')
    r_int = np.arange(0.1,float(sys.argv[2]),0.05)
    
    for i in r_int:
        test = np.interp(i, r, E_rep)
        print(i, test)