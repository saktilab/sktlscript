#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import pymatgen as mg
import pymatgen.util.coord as cu

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')
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
index2print = sorted(index2print)

symbols = []
symbol_list = []
nat_list = []
lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)

#read the first geometry from the trajectory file
with open(traj_file, 'r') as f:
    nat = int(f.readline())
    info = f.readline()
    arr = info.split()
    t1 = float(arr[3])
    for i in range(0, nat):
        line = f.readline()
        arr = line.split()
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
    nat = int(f.readline())
    info = f.readline()
    arr = info.split()
    t2 = float(arr[3])
    dt = t2-t1

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

first = True
first_coords = []

step = 1
fout = open('com.out', 'w')


times = []
msds = []
msds_x = []
msds_yz = []

with open(traj_file, 'r') as f:
    for line in f:            
        summ = 0.0
        summ_x = 0.0
        summ_yz = 0.0
        coords = []
        nat = int(line)
        info = next(f)
        arr = info.split()
print (nat, file=fout)
print ('#Time(ps) msd msd(x) msd(yz)', file=fout)    
        for i in range(0, nat):
            line = next(f)
            if step >= opt.start:
                if (i in index2print):
                    arr = line.split()
                    coords.append( np.array(list(map(float, arr[1:4])),dtype=np.double))
            
        if step < opt.start:
            step += 1
            continue
        if (first):
            first_coords = coords.copy()
            first = False
            previous_coords = coords.copy()
        else:
            for i in range(len(coords)):
                pre_site = previous_coords[i]
                cur_site = coords[i]
                frac_pre = box.get_fractional_coords(pre_site)
                frac_cur = box.get_fractional_coords(cur_site)
                dv = cu.pbc_shortest_vectors(box, frac_pre, frac_cur)[0][0]
                new_site = pre_site + dv
                coords[i] = new_site
                ref_site = first_coords[i]
                dv = new_site - ref_site
                lx = dv[0]*dv[0]
                lyz = (dv[1]*dv[1]+dv[2]*dv[2]) 
                summ_x += lx
                summ_yz += lyz 
                summ += (lx+lyz) 
            msd = (summ/len(first_coords))
            msd_x = (summ_x/len(first_coords))
            msd_yz = (summ_yz/len(first_coords))
            time = dt*step
            step += 1
            print('{:10.4f} {:20.12f} {:20.12f} {:20.12f}'.format(time/1000.0, msd, msd_x, msd_yz), file=fout)
            times.append(time/1000.0)
            msds.append(msd)
            msds_x.append(msd_x)
            msds_yz.append(msd_yz)
            previous_coords = coords



