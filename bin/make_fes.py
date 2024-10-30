#!/usr/bin/env python3
import sys
from scipy.constants import k
import numpy as np
from scipy import integrate
from math import log, exp
s = [] # Collective variables
u = [] # Gassian bias potential
T = float(sys.argv[2])
t = []
F = []

with open(sys.argv[1], 'r') as f:
    for line in f:
        if ('Coordinate' in line):
            arr=line.split()
            s.append(float(arr[2]))
        if ('height' in line):
            arr=line.split()
            u.append(float(arr[3]))
        if ('width' in line):
            arr=line.split()
            width = float(arr[3])
        if ('GAUSSIAN BIAS POTENTIAL' in line):
            arr=line.split()
            t.append(float(sys.argv[3])*int(arr[3])/1000) # time in ps

t = np.array(t,dtype='float',)            
s = np.array(s,dtype='float')
u = np.array(u,dtype='float')
V = np.column_stack((t,u))
np.savetxt('time_bias.dat', V)


F = integrate.cumtrapz(u,t, initial=0)
fes = np.column_stack((s,F))
