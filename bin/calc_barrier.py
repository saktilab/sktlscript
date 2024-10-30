#!/usr/bin/env python3

import numpy as np
from scipy.signal import argrelextrema
import sys
from math import log
kJ2eV = 1.0364e-2   
f = open(sys.argv[1], "r")
dt = np.dtype('float')
x = []
y = []
f.readline()
for line in f:
    arr = line.split()
    x.append(arr[0])
    y.append(arr[1])


x = np.array(x, dtype=dt)
y = np.array(y, dtype=dt)
    
sortid = np.argsort(x)
x = x[sortid]
y = y[sortid]

maxm = argrelextrema(y, np.greater)
minm = argrelextrema(y, np.less)
maxval = [float(i) for i in y[maxm]]
minval = [float(i) for i in y[minm]]

print('Minimum absis: {}'.format(x[minm]))
print('Maximum absis: {}'.format(x[maxm]))
print('Minimum ordinate: {}'.format(y[minm]))
print('Maximum ordinate: {}'.format(y[maxm]))



###Activation barrier for forward and backward reactions

if (len(minval) == 3):
    # If exists two transition states
    Ef1 = (maxval[1]-minval[2])*2625.5 # in kJ/mol
    Ef2 = (maxval[0]-minval[1])*2625.5
    Eb1 = (maxval[0]-minval[0])*2625.5
    Eb2 = (maxval[1]-minval[1])*2625.5
    ##Now calculating delta G
    deltaG = (minval[0]-minval[-1])*2625.5 # in kJ/mol
    print('deltaG: {:.2f} kJ/mol ({:.2f} eV)'.format(deltaG, deltaG*kJ2eV))
    print('deltaG*(backward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef1, Ef1*kJ2eV))
    print('deltaG*(backward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef2, Ef2*kJ2eV))
    print('deltaG*(forward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb1, Eb1*kJ2eV))
    print('deltaG*(forward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb2, Eb2*kJ2eV))
elif (len(minval) == 1):
    Ef1 = (maxval[0]-minval[0])*2625.5 
    Eb1 = (maxval[0]-y[0])*2625.5
    print('deltaG*(forward): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef1, Ef1*kJ2eV))    
    print('deltaG*(backward): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb1, Eb1*kJ2eV))
    deltaG = (y[0]-minval[0])*2625.5
    print('deltaG: {:.2f} kJ/mol ({:.2f} eV)'.format(deltaG, deltaG*kJ2eV))
else:
    Ef1 = (maxval[0]-minval[1])*2625.5 # in kJ/mol
    Ef2 = 0
    Eb1 = (maxval[0]-minval[0])*2625.5
    Eb2 = 0
    ##Now calculating delta G
    deltaG = (minval[0]-minval[-1])*2625.5 # in kJ/mol
    print('deltaG: {:.2f} kJ/mol ({:.2f} eV)'.format(deltaG, deltaG*kJ2eV))
    print('deltaG*(backward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef1, Ef1*kJ2eV))
    print('deltaG*(backward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Ef2, Ef2*kJ2eV))
    print('deltaG*(forward1): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb1, Eb1*kJ2eV))
    print('deltaG*(forward2): {:.2f} kJ/mol ({:.2f} eV)'.format(Eb2, Eb2*kJ2eV))
# Eb = (maxval[0]-minval[0])*2625.5 # in kJ/mol
# print('Forward barrier: {:.2f} kJ/mol'.format(Ef))
# print('Backward barrier: {:.2f} kJ/mol'.format(Eb))

