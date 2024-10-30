#!/usr/bin/env python3

import sys

Energy_LV = []
Energy_LP = []

with open("orca.log", 'r') as f:
    for line in f:
        if ("LV ( 1) V  1") in line:
            arr=line.split()
            if len(arr) > 13:
                Energy_LV.append(float(arr[-3]))
        if ("LV ( 2) V  1") in line:
            arr=line.split()
            if len(arr) > 13:
                Energy_LV.append(float(arr[-3]))
        if ("LV ( 3) V  1") in line:
            arr=line.split()
            if len(arr) > 13:
                Energy_LV.append(float(arr[-3]))
        if ("LP ( 1) V  1") in line:
            arr=line.split()
            if len(arr) > 14:
                Energy_LP.append(float(arr[-3]))
print("Lone Vacant Energy = {:.2f} kcal/mol".format(sum(Energy_LV)))
print("LV_Energy/Bond = {:.2f} kcal/mol".format(sum(Energy_LV)/len(Energy_LV)))
print("Lone Pair Energy = {:.2f} kcal/mol".format(sum(Energy_LP)))
print("LP_Energy/Bond = {:.2f} kcal/mol".format(sum(Energy_LP)/len(Energy_LP)))



