#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pymatgen as mg
import copy
import statistics as st

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('nac', metavar='nac', type=str)
parser.add_argument('-l', '--latt', type=str, required=True)
# parser.add_argument('-c', '--charge', type=float, default=-0.60)
parser.add_argument('-s', '--stepend', type=int, default=100)
parser.add_argument('-nmin', type=int, default=1, help='number of minimum charge to list')
parser.add_argument('-nmax', type=int, default=1, help='number of minimum charge to list')


group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-e', '--elem', type=str, nargs='+', help='Elements to be analyzed')
group.add_argument('-i', '--index', nargs='+', type=int, help='Index to be analyzed')
group.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')
opt = parser.parse_args(sys.argv[1:])


traj_file = opt.traject
nac_file = opt.nac

lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)


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


#read the first geometry from the trajectory file
with open(traj_file, 'r') as f:
    nat = int(f.readline())
    info = f.readline()
    arr = info.split()
    t1 = float(arr[3])
    # t1 = 0 
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


#parse the mulliken charge file and compute the net atomic charge(NAC)
step = 1
max_charge_atoms = []
min_charge_atoms = []
charge_thresh = []
# avg_charge = []
with open(nac_file, 'r') as f:
    for line in f:
        nat = int(line)
        stepline = next(f)
        arr = stepline.split()
        step += 1
      
        charge_index = []
        charge = []
        for ind, elem in enumerate(symbols):
            line = next(f)
            if (ind in index2print):
                arr = line.split()
                #nac = float(arr[1])
                nac=float(arr[2])
                charge_index.append([ind, nac])
                charge.append(nac)
        sorted_charge_index = sorted(charge_index, key=lambda x: x[1])
        max_index = len(sorted_charge_index)-opt.nmax-1
        
        max_charges = sorted_charge_index[len(sorted_charge_index)-1:max_index:-1]
        min_charges = sorted_charge_index[0:opt.nmin]
        avg_charge = st.mean(charge)
        charge = np.array(charge, dtype='float')
        sigma = np.std(charge)
        # Use dynamical threshold to determine the proton transfer events
        thresh = avg_charge+5*sigma
        
         
        filtered_max_charge = []
        for item in max_charges:
            if (item[1] < thresh):
            # if (item[1] < hydronium_charge):
                filtered_max_charge.append([-1, item[1]])
            else:
                filtered_max_charge.append(item)
        max_charge_atoms.append(filtered_max_charge)
        min_charge_atoms.append(min_charges)
        # avg_charge.append(sum([x[0] for x in charge_index[1:-1]])/float(len(charge_index)-2))
        
  

fout = open('charge_list.dat', 'w')
header_format_str = ' {:>10s} {:>10s}'*(opt.nmin+opt.nmax)
header_text = ['MaxIndex', 'MaxCharge']*(opt.nmax) + ['MinIndex', 'MinCharge']*(opt.nmin)
row_format_str = ' {:>10d} {:10.4f}'*(opt.nmin+opt.nmax)
print (('#{:>10s}' + header_format_str).format('Time(ps)', *header_text))
print (('#{:>10s}' + header_format_str).format('Time(ps)', *header_text), file=fout)


flatten = lambda l: [item for sublist in l for item in sublist]


for i in range(0, len(max_charge_atoms)):
    time = dt*i
    print (('{:06d} {:10.4f}'+row_format_str).format(i, time/1000.0, *flatten(max_charge_atoms[i]), *flatten(min_charge_atoms[i])  ), file=fout)
    print (('{:06d} {:10.4f}'+row_format_str).format(i, time/1000.0, *flatten(max_charge_atoms[i]), *flatten(min_charge_atoms[i])  ))    

print('*'*80)



            

