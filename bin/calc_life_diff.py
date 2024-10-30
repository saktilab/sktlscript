#!/usr/bin/env python3

import sys
from scipy.constants import k, h, R
from math import exp
f = open(sys.argv[1], 'r')
T = float(sys.argv[2])
l = float(sys.argv[3])
kappa = 1.0  
pmf = []
# Extract the proton transfer barrier
for line in f:
    arr=line.split()
    pmf.append(arr[1])
Ea = float(pmf[0])

#Calculate the life time
rate = kappa*(k*T/h)*exp(-Ea*4.184*1000/(R*T))  
tau = 1.0/rate*1e12 # in picoseconds
D = l**2/(2*tau) # in Angstrom^2/ps
print("deltaF* = {:.4} kcal/mol".format(Ea))
print("Life time(tau) = {:.4} ps".format(tau))
print("D = {:.4} Angstrom^2/ps".format(D))
if (D<1.0e-6):
    print("Need metadynamics simulation to confirm the accuracy!")