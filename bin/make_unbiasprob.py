#!/usr/bin/env python3

import sys
import numpy as np
from scipy.constants import k

kT = k*2.294e+17*float(sys.argv[2]) # in hartree

F = []
cv = []
with open(sys.argv[1]) as f:
    for line in f:
        arr = line.split()
        cv.append(float(arr[0]))
        F.append(float(arr[1]))

cv = np.array(cv, dtype='float')
F = np.array(F, dtype='float')
P = np.exp(-F/kT)
unbiasprob = np.column_stack((cv,P))
np.savetxt("unbiasprob.dat", unbiasprob)