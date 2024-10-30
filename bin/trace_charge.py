#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pymatgen as mg
import copy

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('nac', metavar='nac', type=str)
parser.add_argument('-l', '--latt', type=str, required=True)
parser.add_argument('-c', '--charge', type=float, default=-0.65)
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
fout1 = open("number_of_species.dat", "w")
fout2 = open("Index_of_oh.dat", "w")
fout3 = open("transfer_event.dat", "w")
print("{:2s} {:2s} {:2s} {:2s}".format('# Time [ps]', 'N(OH-)' ,'N(water)', 'N(PT events)'),file=fout1)
# print("{:2s} {:2s}*len(index_oh)".format('# Time [ps]', "Index_oh"),file=fout2)
# avg_charge = []
with open(nac_file, 'r') as f:
    count_oh = []
    for line in f:
        nat = int(line)
        stepline = next(f)
        arr = stepline.split()
        step += 1

        charge_index = []
        serial_oh = []
        serial_wat = []
        
        for ind, elem in enumerate(symbols):
            line = next(f)
            if (ind in index2print):
                arr = line.split()
                # nac = float(arr[1])
                nac=float(arr[2])
                charge_index.append([ind, nac])
        
        for ind, nac in charge_index:
            if (nac < -0.75):
                serial_oh.append(ind)
            else:
                serial_wat.append(ind)
        
        print("{:2f} {:2d} {:2d}".format(step/200.0,len(serial_oh),len(serial_wat)),file=fout1)
        count_oh.append(len(serial_oh))
    PT = [abs(j-i) for i,j in zip(count_oh[:-1], count_oh[1:])]
    for i in PT:
        print(i, file=fout3)
        # print(("{:2f}" ("{}"*len(serial_oh))).format(step/200.0,*serial_oh),file=fout2)

# Now we identify the neighbor of each index

# final_pairs = []

# print ()
# print ('Pairs to be printed:')
# for li in list_pairs:
    
#     if (li[0] < len(symbols) and li[1] < len(symbols)):
#         print ('{:5d} {:5d}'.format(*li))
#         final_pairs.append([li[0], li[1]])
# print ()
# step = 0
# format_str = ' {:>15s}'*len(final_pairs)
# headers = ['{}({})-{}({})'.format(symbols[i], i, symbols[j], j) for i, j in final_pairs ]
# print(('{:>10s} ' + format_str + ' {:15s}').format('#Time(ps)', *headers, 'Average'), file=fout)
pbc = False
if (len(lattice) == 3):
    pbc = True

# with open(traj_file, 'r') as f:
#     for line in f:
#         time = dt*step
#         step+=1
#         nat = int(line)
#         next(f)
#         coords = []
#         for i in range(0, nat):
#             line = next(f)
#             arr = line.split()
#             coords.append(list(map(float, list(arr[1:4]))))
#         if (pbc):
#             my_str = mg.Structure(box, species=symbols, coords=coords, coords_are_cartesian=True)
#         else:
#             my_str = mg.Molecule(species=symbols, coords=coords)
        # for i in serial_oh:
        #     print(i)
        # # my_str.get_neighbors()
            # neighbors = my_str.get_neighbors(site, 2.5)
            # print(neighbors)
         # distances = []
        # for ind_1, ind_2 in final_pairs:
        #     distance = my_str.get_distance(ind_1, ind_2)
        #     distances.append(distance)
        # format_str = ' {:15.4f}'*len(distances)
        # print(('{:>10.4f} ' + format_str + ' {:15.4f}').format(time, *distances, sum(distances)/len(distances)),file=fout)
        # print(('{:>10.4f} ' + format_str + ' {:15.4f}').format(time, *distances),file=fout)    

