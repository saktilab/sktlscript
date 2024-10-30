#!/usr/bin/env python3 
import sys
from scipy.constants import k, pi

epsilon = 2.837297
eta = 0.896e-3
f = open(sys.argv[1], 'r')
for line in f:
    if  "Diffusion coefficient (1/6)" in line:
        arr = line.split()
        D_pbc = float(arr[3])*1e-8
        break

with open(sys.argv[2], 'r') as f:
    for line in f:
        if 'TV' in line:
            arr = line.split()
            L = float(arr[1])*1e-10
            break
T = float(sys.argv[3])
corr =  k*T*epsilon/(6*pi*eta*L)
D0 = (D_pbc + corr)*1e8
print(D0)