#!/usr/bin/env python3


import sys
x = []
y = []
z = []
with open(sys.argv[1],'r') as f:
	for line in f:
		if "POINTS" in line:
			for i in range(3):
				arr = next(f).split()
				x.append(arr[0])
				y.append(arr[1])
				z.append(arr[2])
x_frame = []
y_frame = []
z_frame = []
Natom = 0
with open(sys.argv[2],'r') as f:
	for line in f:
		if "POINTS" in line:
			for i in range(Natom):
				arr = next(f).split()
				x_frame.append(arr[0])
				y_frame.append(arr[1])
				z_frame.append(arr[2])

ads_symb = []
Natom = 0
with open(sys.argv[3],'r') as f:
	for line in f:
		if "Number of atoms" in line:
			Natom += int(next(f)) 
		if "atomic positions" in line:
			for i in range():
				arr = next(f).split()
				ads_symb.append(arr[1][0])
Nat_frame = 0
frame_symb = []
with open(sys.argv[4],'r') as f:
	next(f)
	next(f)
	next(f)
	Nat_frame += int(next(f))
	for i in range(Nat_frame):
		arr = next(f).split()
		frame_symb.append(arr[4])
	

symb = frame_symb + ads_symb
x = x + x_frame
y = y + y_frame
z = z + z_frame
print(x_frame)
# for symb,x,y,z in zip(symb,x,y,z):
# 	print(symb,x,y,z)



