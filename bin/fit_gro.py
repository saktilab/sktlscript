#!/usr/bin/env python3
import numpy as np
from scipy.optimize import leastsq

file_gro = open("gro.dat", "r") 
h = []
time = []
next(file_gro)
for line in file_gro:
    arr = line.split()
    time.append(arr[0])
    h.append(float(arr[1]))
time = np.array(time, dtype=float)
h = np.array(h, dtype=float)
def lin_fit(x, y):
    ''' 
    Fits a linear fit of the form y=ax
    '''
    fitfunc = lambda params, x: params[0] * x
    errfunc = lambda p, x, y: fitfunc(p,x) -y

    init_a = 0.5
    init_p = np.array(init_a)

    p1, success = leastsq(errfunc, init_p.copy(), args=(x,y))
    f = fitfunc(p1, x)
    return p1, f

fitting = lin_fit(time,h)
diff = float(fitting[0]*2.5**2/6)
print("Diffusion const.: {:2f} squared-angstrom/ps".format(diff))
print("PT rate: {:2f} /ps".format(float(fitting[0])))
