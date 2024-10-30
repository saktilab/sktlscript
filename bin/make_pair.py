#!/usr/bin/env python3

import sys

file_ind_H = open(sys.argv[1], "r")
file_ind_O = open(sys.argv[2], "r")

indexes_H = []
for line in file_ind_H:
    indexes_H.append(int(line))
indexes_O = []
for line in file_ind_O:
    indexes_O.append(int(line))

# fout = open("pair_H_O.dat","w")
for i in indexes_H:
    fout = open("pair_H_O_%i.dat" %i, "w")
    for ind_O in indexes_O:
        print(i,ind_O, file=fout)
