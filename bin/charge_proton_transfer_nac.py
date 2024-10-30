#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('nac', metavar='nac', type=str)
parser.add_argument('-p', '--proton', type=bool, default=False,  help='Protonated system?')
parser.add_argument('-c', '--charge', type=float, default=-0.73, help='charge limit for hydronium')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-e', '--elem', type=str, nargs='+', help='Elements to be analyzed')
group.add_argument('-i', '--index', nargs='+', type=int, help='Index to be analyzed')
group.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')
opt = parser.parse_args(sys.argv[1:])


traj_file = opt.traject
nac_file = opt.nac


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
avg_charge = []
with open(nac_file, 'r') as f:
    for line in f:
        nat = int(line)
        stepline = next(f)
        arr = stepline.split()
        step += 1
        
        charge_index = []
        for ind, elem in enumerate(symbols):
            line = next(f)
            if (ind in index2print):
                arr = line.split()
                nac = float(arr[2])
                charge_index.append([nac, ind])
                
        sorted_charge_index = sorted(charge_index)
        if (opt.proton) :
            if (sorted_charge_index[-1][0] < opt.charge):
                prev_index = max_charge_atoms[-1][1]
                prev_array_index = 0
                for i, item in enumerate(charge_index):
                    if (charge_index[1] == prev_index):
                        prev_array_index = i
                        break
                max_charge_atoms.append([charge_index[prev_array_index][0], prev_index])
                print('Warning: hydronium gone!?')
            else:
                max_charge_atoms.append(sorted_charge_index[-1])
        else:
            max_charge_atoms.append(sorted_charge_index[-1])
        min_charge_atoms.append(sorted_charge_index[0])
        avg_charge.append(sum([x[0] for x in charge_index[1:-1]])/float(len(charge_index)-2))
        

num_pros = 0

iList = max_charge_atoms[0][1]
jList = max_charge_atoms[0][1]
if (max_charge_atoms[1][1] != max_charge_atoms[0][1]):
    num_pros = 1
    iList = max_charge_atoms[1][1]

lList = max_charge_atoms[0][1]
time = dt*1
print ('#{:>9s} {:>8s} {:>8s} {:>10s} {:>8s} {:>10s} {:>10s}'.format('Time(ps)', '#Proton', 'MaxIndex', 'MaxCharge', 'MinIndex', 'MinCharge', 'AvgCharge'))
print ('{:10.4f} {:>8d} {:>8d} {:>10.4f} {:>8d} {:>10.4f} {:>10.4f}'.format(time/1000.0, num_pros, max_charge_atoms[1][1], max_charge_atoms[1][0], min_charge_atoms[1][1], min_charge_atoms[1][0], avg_charge[1]   ))
for i in range(2, len(max_charge_atoms)):
    time = dt*i
    if (max_charge_atoms[i][1] == max_charge_atoms[i-1][1]):
        pass
    else:
        if (max_charge_atoms[i][1] != iList and max_charge_atoms[i][1] == jList):
            if (iList == lList):
                num_pros += 1
            else:
                num_pros -= 1
                lList = max_charge_atoms[i][1]
        else:
            num_pros += 1
        jList = iList
        iList = max_charge_atoms[i][1]

    print ('{:10.4f} {:>8d} {:>8d} {:>10.4f} {:>8d} {:>10.4f} {:>10.4f}'.format(time/1000.0, num_pros, max_charge_atoms[i][1], max_charge_atoms[i][0], min_charge_atoms[i][1], min_charge_atoms[i][0], avg_charge[i]  ))

with open('id_hydronium.dat', 'w') as f:
    for charge, atom in max_charge_atoms:
        print ('serial {}'.format(atom+1), file=f)

with open('id_hydroxide.dat', 'w') as f:
    for charge, atom in min_charge_atoms:
        print ('serial {}'.format(atom+1), file=f)
    

