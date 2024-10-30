#!/usr/bin/env python3
import sys
import argparse
import decimal
import numpy as np
    


def drange(x, y, jump):
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)

d = """
Calculating electric double layer, namely, the voltage as a function of the distance of 
electrode-electrolyte interface
"""

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
parser.add_argument("-geo", type=str, help="Geometry file")
parser.add_argument("-cha", type=str, help="Charge file")

args = parser.parse_args()

element = 1.60217657
coulomb = 8.9876
x = []
y = []
z = []
symbol = []
charge = []
# Obtaining coordinate information
with open(args.geo, "r") as f:
    natom = int(next(f))
    info = next(f)
    for line in f:
        arr = line.split()
        symbol.append(arr[0])
        x.append(arr[1])
        y.append(arr[2])
        z.append(arr[3])

x = np.array(x,dtype='float')
y = np.array(y,dtype='float')
z = np.array(z,dtype='float')

#Obtaining charge information
with open(args.cha, "r") as f:
    next(f)
    next(f)
    for line in f:
        arr = line.split()
        charge.append(arr[1]) 

charge = np.array(charge, dtype='float')

#Calculate distance
xmax = np.amax(x)+1
xmin = np.amin(x)-1
ymax = np.amax(y)+1
ymin = np.amin(y)-1
zmax = np.amax(z)+1
zmin = np.amin(z)-1

delta_x = (xmax - xmin)*10
delta_y = (ymax - ymin)*10
delta_z = (zmax - zmin)*10

potential = 0
potential_x = 0
for i in drange(1, delta_x+1,1):
    rx = 0
    rx = xmin + (i - 1) * 0.1
    potential_y = 0
    for j in drange(1,delta_y+1,1):
        ry = 0
        ry = ymin + (j-1) * 0.1
        potential_z = 0
        for k in drange(1, delta_z+1, 1):
            rz = 0
            rz = zmin + (k-1) * 0.1
            for l in range(1, natom):
                distance = 0
                pot = 0
                distance = ((rx - x[l])**2 + (ry - y[l])**2 + (rz - z[l])**2)**0.5
                # potential = 0
                if (distance != 0):
                    pot = element * coulomb * charge[l] / distance
                    potential = potential + pot
                if (l == natom):
                    potential_z[k] = potential
                    potential = 0
            if (k == delta_z + 1):
                potential_y[j] = np.sum(potential_z)
        if (j == delta_y + 1):
            potential_x= np.sum(potential_y)/((potential_y+1) * (potential_z+1))
    print(i, potential_x) 
 
