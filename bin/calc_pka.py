#!/usr/bin/env python3

import numpy as np
import os
from scipy.signal import argrelextrema
import sys
from math import log

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

##Now calculating delta G
#temp = os.system("~/bin/md_stats.py -p temp dftb.out | tail -n 2 | head -n 1 |awk '{print $2}'")
#print(temp)
deltaG = (minval[0]-minval[-1])*2625.5 # in kJ/mol
#pKa = deltaG*1000/(2.303*8.314*float(temp))
pKa = deltaG*1000/(2.303*8.314*float(sys.argv[2]))
# print('{:.2f}'.format(pKa))

###Activation barrier for forward and backward reactions
Ef = (maxval[0]-minval[-1])*2625.5 # in kJ/mol
Eb = (maxval[0]-minval[0])*2625.5 # in kJ/mol
# print('Forward barrier: {:.2f} kJ/mol'.format(Ef))
# print('Backward barrier: {:.2f} kJ/mol'.format(Eb))
print('deltaG: {:.2f} kJ/mol'.format(deltaG))
print('deltaG*(forward): {:.2f} kJ/mol'.format(Ef))
print('deltaG*(backward): {:.2f} kJ/mol'.format(Eb))
print('pKa: {:.2f} '.format(pKa))
