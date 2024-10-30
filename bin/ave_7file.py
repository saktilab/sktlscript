#!/usr/bin/env python3
import sys
import numpy as np
import os
f1 = np.loadtxt(sys.argv[1], dtype='float')
f2 = np.loadtxt(sys.argv[2], dtype='float')
f3 = np.loadtxt(sys.argv[3], dtype='float')
f4 = np.loadtxt(sys.argv[4], dtype='float')
f5 = np.loadtxt(sys.argv[5], dtype='float')
f6 = np.loadtxt(sys.argv[6], dtype='float')
f7 = np.loadtxt(sys.argv[7], dtype='float')
f_ave = (f1 + f2 + f3 + f4 + f5 + f6 + f7)/7.0
# print(f_ave) 
# f_ave.tofile("average.dat")
np.savetxt("average.dat", f_ave)