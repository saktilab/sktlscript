#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='List the charge during trajectory')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('mulliken', metavar='mulliken', type=str)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-e', '--elem', type=str, nargs='+', help='Elements to be analyzed')
group.add_argument('-i', '--index', nargs='+', type=int, help='Index to be analyzed')
group.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')
opt = parser.parse_args(sys.argv[1:])


#n_orbs : number of shells. 
n_orbs = {
    'Pt': 3,
    'Br': 3,
    'O': 2,
    'C': 2,
    'N': 2,
    'H': 1,
    'S': 3
}

net_charge = {
    'Pt': 10,
    'Br': 7,
    'O': 6,
    'C': 4,
    'N': 5,
    'H': 1,
    'S': 6
}


traj_file = opt.traject
mull_file = opt.mulliken


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
all_charge_index = []
with open(mull_file, 'r') as f:
    for line in f:
        nat = int(line)
        stepline = next(f)
        arr = stepline.split()
        # print ('STEP {}'.format(step))
        step += 1
        
        charge_index = []
        
        for ind, elem in enumerate(symbols):
            norb = n_orbs[elem]
            total_pop = 0.0
            for j in range(0, norb):
                line = next(f)
                arr = line.split()                
                total_pop += float(arr[1])
            if (ind in index2print):
                nac = net_charge[elem]-total_pop
                charge_index.append([ind, nac])
        
        all_charge_index.append(charge_index)           
        


indexes = [x[0] for x in all_charge_index[0]]
num_indexes = len(indexes)


fout = open('charge_list.out', 'w')

formatstr_header = '#{:>9s}' + ' {:>10d}'*num_indexes
print (formatstr_header.format('Time(ps)', *indexes), file=fout)
step = 0
formatstr_line = ' {:>10.4f}'*num_indexes
for step_index in all_charge_index:
    values = [x[1] for x in step_index]
    time = dt*step
    step+=1
    line = formatstr_line.format(*values)
    print ('{:10.4f} {}'.format(time/1000.0, line), file=fout)




