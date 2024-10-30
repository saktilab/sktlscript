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
# avg_charge = []
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
                #nac = float(arr[1])
                nac=float(arr[2])
                charge_index.append([ind, nac])
                
        sorted_charge_index = sorted(charge_index, key=lambda x: x[1])
        max_index = len(sorted_charge_index)-opt.nmax-1
        
        max_charges = sorted_charge_index[len(sorted_charge_index)-1:max_index:-1]
        min_charges = sorted_charge_index[0:opt.nmin]
        filtered_max_charge = []
        for item in max_charges:
            if (item[1] < opt.charge):
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
    time = dt*i+t1
    print (('{:06d} {:10.4f}'+row_format_str).format(i, time/1000.0, *flatten(max_charge_atoms[i]), *flatten(min_charge_atoms[i])  ), file=fout)
    print (('{:06d} {:10.4f}'+row_format_str).format(i, time/1000.0, *flatten(max_charge_atoms[i]), *flatten(min_charge_atoms[i])  ))    


print('*'*80)


# fouts = []
# for i in range(len(max_charge_atoms[0])):
#     fouts.append(open('id_max_charge_{}.dat'.format(i+1), 'w'))
# for line_step in max_charge_atoms:
#     line = [x[0] for x in line_step]
#     for i in range(len(line)):
#         if  (line[i] == -1):
#             print ("none", file=fouts[i])
#         else:
#             print('index {}'.format(line[i]), file=fouts[i])
# for ff in fouts:
#     ff.close()

# index_max_charges = []


# for i in range(0, len(max_charge_atoms)):
#     index_max_charges.append([x for x in max_charge_atoms[i]])


# pre_max_coords = []
# step = 0
# final_index_max_charges = []
# with open(traj_file, 'r') as f:
#     for line in f:
#         nat = int(line)
#         info = next(f)
#         arr = info.split()
     
#         cur_max_coords = [None]*len(index_max_charges[step])
#         for i in range(0, nat):
#             line = next(f)
#             arr = line.split()

#             if (i in [x[0] for x in index_max_charges[step]]):
#                 for k in range(len(index_max_charges[step])):                   
#                     if (i == index_max_charges[step][k][0]):
#                         cur_max_coords[k] =list(map(float, arr[1:4])) 
#                         break

#         print (cur_max_coords, index_max_charges[step])
#         #sort index according to distance
#         if (step ==0):
#             pre_max_coords = cur_max_coords
#             final_index_max_charges.append(index_max_charges[step])
#         else:
#             sorted_index_max = []
            
#             temp_max_map = {}
#             for i, v in zip(index_max_charges[step], cur_max_coords):
#                 if (i[0] > -1):
#                     temp_max_map[i[0]] = [v, i[1]]
            
#             if (len(temp_max_map) == 0):
#                 final_index_max_charges.append(index_max_charges[step])

            
#             final_max_coords = []
#             # print (pre_max_coords, index_max_charges[step])
#             updated_pre_max_coords = []

#             print(pre_max_coords)
#             for i in range(len(pre_max_coords)):
#                 if (pre_max_coords[i] is None):
#                     continue
#                 # print(temp_max_map)
#  #               print('Find dis of {}'.format(index_max_charges[step-1][i]))
#                 fcoord1 = box.get_fractional_coords(pre_max_coords[i])           
#                 dis_index = []
#                 for k in temp_max_map.keys():
#                     fcoord2 = box.get_fractional_coords(temp_max_map[k][0])            
#                     dis, image = box.get_distance_and_image(fcoord1, fcoord2)
#                     dis_index.append([dis, k, temp_max_map[k][1]])
#  #                   print(' - {}: r: {}, c={}'.format(k, dis, temp_max_map[k][1]))
                
#                 sorted_dis_index = sorted(dis_index)
#                 print(sorted_dis_index)
                
#                 final_index = sorted_dis_index[0][1]                
#                 sorted_index_max.append([final_index, sorted_dis_index[0][2]])
#                 final_max_coords.append(temp_max_map[final_index][0])
#                 updated_pre_max_coords.append(temp_max_map[final_index][0])                
#                 del temp_max_map[final_index]
#                 # print(temp_max_map)

            
#             final_index_max_charges.append(sorted_index_max)
# #            print([x[0] for x  in sorted_index_max])
#             pre_max_coords = updated_pre_max_coords


#         step+=1

# fouts = []
# for i in range(len(final_index_max_charges[0])):
#     fouts.append(open('sorted_charge_list_{}.dat'.format(i+1), 'w'))
# for step, line_step in enumerate(final_index_max_charges):
#     time = (dt*step+t1)/1000.0
#     for i in range(len(line_step)):
#         print('{:10.4f} {:5d} {:10.4f}'.format(time, *line_step[i]), file=fouts[i])
        
# for ff in fouts:
#     ff.close()
            

