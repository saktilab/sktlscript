#!/usr/bin/env python3

f = open('biaspot', 'r')
g = open('colvar.dat', 'r')

for line in f:
    gaus_step = line.split()
    arr = line.split()
    print(gaus_step)
