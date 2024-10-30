#!/usr/bin/env python3

import sys  
import argparse
import numpy as np
import pymatgen as mg
import matplotlib
matplotlib.use('PS')
from sklearn import linear_model
from matplotlib import pyplot as plt
import pymatgen.util.coord_utils as cu
import math
from scipy import integrate


parser = argparse.ArgumentParser(description='Velocity autocorelation function')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-s', '--start', type=int, default=0, help='StartTime')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-e', '--elem', type=str, nargs='+', help='Elements to be analyzed')
group.add_argument('-i', '--index', nargs='+', type=int, help='Index to be analyzed')
group.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')
opt = parser.parse_args(sys.argv[1:])


traj_file = opt.traject
print(opt.start)

index2print = []

if (opt.elem):
    print()
    print('Elements to be analyzed: ', ", ".join(opt.elem))
elif (opt.index):
    index2print = opt.index
elif (opt.index_file):
    with open(opt.index_file, 'r') as f:
        for line in f:
            line_list = list(map(int, line.split()))
            index2print += line_list

symbols = []
symbol_list = []
nat_list = []

veloc_all_step = []

# fout2 = open('oveloc.dat', 'w')
#read the first geometry from the trajectory file
step = 0
t1 = 0.0
t2 = 0.0
delta_t = 0.0
with open(traj_file, 'r') as f:
    for line in f:
        veloc = []
        nat = int(line)
        info = next(f)
        if (step == 0):
            arr = info.split()
            t1 = float(arr[3])
        elif (step == 1):
            arr = info.split()
            t2 = float(arr[3])
            delta_t = t2-t1
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            if (step == 0):
                symbols.append(arr[0])           
                if (opt.elem and (arr[0] in opt.elem )):
                    index2print.append(i)
                if (len(symbol_list) == 0):
                    symbol_list.append(arr[0])
                    nat_list.append(1)
                else:
                    if (arr[0] == symbol_list[-1]):
                        nat_list[-1] += 1
                    else:
                        symbol_list.append(arr[0])
                        nat_list.append(1)
            if (i in index2print):
                veloc.append(np.array(list(map(float, arr[1:4])),dtype=np.double))
                # print('{} {} {} {}'.format(*arr[0:4]),file=fout2)
        veloc_all_step.append(veloc)  
        step += 1
    

index_list = [index2print[x:x+5] for x in range(0, len(index2print),5)]


print ()
print ('Indexes to be analyzed:')
num = 0
for li in index_list:
    print ('{:3d}:'.format(num*5+1), ('  {:5d}'*len(li)).format(*li))
    num+=1
print ()

if (not opt.index_file):
    with open('index_auto.dat', 'w') as f:
        for li in index_list:
            print (('  {:5d}'*len(li)).format(*li), file=f)

fout = open('vacf.out', 'w')
print ('#Time(ps) z(t)', file=fout)






nsteps = len(veloc_all_step)
vacf = np.zeros(nsteps, dtype=np.double) 
natoms = len(veloc_all_step[0])

print(nsteps, natoms)
for i in range(0, nsteps):
    vac_i = 0.0
    for j in range(natoms):
        norm1 = sum(veloc_all_step[0][j]**2)**0.5
        norm2 = sum(veloc_all_step[i][j]**2)**0.5
        if (norm1 > 0.0 and norm2 > 0.0):
            vac_i += sum(veloc_all_step[0][j]*veloc_all_step[i][j])/norm1/norm2
    vacf[i] = (vac_i/float(natoms))
vacf /= vacf[0]

times = []
for t in range(0, nsteps):
    value = vacf[t]
    time = t*delta_t
    times.append(time/1000.0)
    print('{:12.6f} {:20.12f}'.format(time/1000.0, value ), file=fout)

print(integrate.simps(abs(vacf), times)/3.0)

# print ((1/3)*(delta_t/1000.0)*sum(list(map(abs, vacf))))


