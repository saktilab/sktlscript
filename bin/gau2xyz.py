#!/usr/bin/env python3

import sys

symb = []
z = []
rx = []
ry = []
rz = []

with open(sys.argv[1],'r') as f:
	for line in f:
		arr = line.split()
		z.append(arr[0])
		rx.append(arr[2])
		ry.append(arr[3])
		rz.append(arr[4])

fname = sys.argv[1].split(".")

z2zsymb = {"1":"H","8":"O","6":"C","40":"Zr","9":"F","17":"Cl"}

for i in z:
	symb.append(z2zsymb[i])

with open("{}.xyz".format(fname[0]),'w') as fout:
	print("{}".format(len(z)),file=fout)
	print("",file=fout)
	for symb,x,y,z in zip(symb,rx,ry,rz):
		print(symb,x,y,z,file=fout)

