#!/usr/bin/env python3

import sys
import pymatgen as mg
import glob

diss_all = []

elem_1 = sys.argv[1]
elem_2 = sys.argv[2]

max_limit = 12.0

filenames = []
for name in sys.argv[3:]:
    filenames += glob.glob(name)

for name in sys.argv[3:]:
    ala = mg.Structure.from_file(name)
    nat = len(ala.sites)

    diss = []

    for i in range(nat):
        if (str(ala.sites[i].specie) != elem_1):
            continue

#        print('Atom : ', i , ala.sites[i])
        
        lista = ala.get_neighbors(ala.sites[i], 8.0, include_index=True)
        sorted_lista = sorted(lista, key=lambda x: x[1])
        listb = [x for x in sorted_lista if str(x[0].specie) == elem_2]

        if (len(listb) == 0):
            raise ValueError('Pair ', elem_1, '-', elem_2, ' not found')
        min_dis = listb[0][1]
        diss.append(min_dis)
    if (len(diss) > 0):
        sort_dis = sorted(diss)
    #    print(sort_dis)
        diss_all.append(sort_dis[0])
        diss_all.append(sort_dis[-1])

sort_dis_all = sorted(diss_all)
if (len(sort_dis_all) > 0):
    print(sort_dis_all[0]/0.529177, sort_dis_all[-1]/0.529177)
