#!/usr/bin/env python3
import sys
states_1 = []
states_2 = []
energy_1 = []
energy_2 = []
thresh = 0.001
skip = ['KPT']
with open(sys.argv[1], 'r') as f:
    f = filter(lambda x: x.strip(), f)
    for line in f:
        if not any(x in line for x in skip):  
            arr = line.split()
            state = int(arr[0])
            energy = float(arr[1])
            occ = float(arr[2])
            if occ < thresh:
                energy_2.append(energy)
                states_2.append(state)
            if occ > thresh: 
                energy_1.append(energy)
                states_1.append(state)
    state_cbm = states_2[energy_2.index(min(energy_2))]
    energy_cbm = min(energy_2)
    state_vbm = states_1[energy_1.index(max(energy_1))]
    energy_vbm = max(energy_1)
    e_gap = energy_cbm - energy_vbm
    print(state_vbm,state_cbm,e_gap)