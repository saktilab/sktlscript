#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pymatgen as mg
import copy

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('chargelist', metavar='chargelist', type=str)
parser.add_argument('-l', '--latt', type=str, required=True)
parser.add_argument('-i', '--index', type=int, required=True)
opt = parser.parse_args(sys.argv[1:])


traj_file = opt.traject
charge_list = opt.chargelist

lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)


def read_xyz(fin):
    coords = []
    line = next(fin)
    nat = int(line)
    next(fin)
    for i in range(nat):
        line = next(fin)
        arr = line.split()
        coords.append(list(map(float, arr[1:4])))
    return coords



symbols = []
symbol_list = []
nat_list = []

charge_histrory = []

ftrj = open(traj_file)

previous_index = opt.index


with open(charge_list) as f:
    for line in f:
        if (line.startswith('#')):
            continue
        arr = line.split()
        step = int(arr[0])
        time = float(arr[1])

        coords = read_xyz(ftrj)
        previous_coords = coords[previous_index]

        distances = []

        no_compare = False
        fcoord1 = box.get_fractional_coords(previous_coords)
        for i in range((len(arr)-2)//2):
            index = int(arr[2+i*2])
            charge = float(arr[3+i*2])
            
            # no oxygen, then use the previous indexes
            if (index == -1 and i == 0):
                # last_item = charge_histrory[-1].copy()
                last_item = [0,0.0,-1,0.0,[0.0,0.0,0.0]]
                last_item[0] = step
                last_item[1] = time
                last_item[2] = -1
                last_item[3] = 0.0
                last_item[4] = [0.0, 0.0, 0.0]
                charge_histrory.append(last_item)             
                no_compare = True 
                break 
            elif (index == -1):
                break
            else:
                fcoord2 = box.get_fractional_coords(coords[index])  
                dis, image = box.get_distance_and_image(fcoord1, fcoord2)
                distances.append([dis, index, charge, fcoord2])
                
        if (no_compare):
            continue
        else:
            sorted_distances = sorted(distances)
            charge_histrory.append([step, time, sorted_distances[0][1], sorted_distances[0][2], sorted_distances[0][3]])

print ('# {:4s} {:>10s} {:5s} {:>10s} {:>12s} {:>12s} {:>12s}'.format('Step', 'Time(ps)', 'Index', 'Charge', 'Coords(x)', 'Coords(y)' , 'Coords(z)'))
for item in charge_histrory:
    print ('{:6d} {:10.4f} {:5d} {:10.4f} {:12.6f} {:12.6f} {:12.6f}'.format(item[0], item[1], item[2], item[3], *item[4]))
            

ftrj.close

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
            

