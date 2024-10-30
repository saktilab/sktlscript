#!/usr/bin/env python3
import sys
import numpy as np
import os
# f1 = np.loadtxt(sys.argv[1], dtype='float')
# f2 = np.loadtxt(sys.argv[2], dtype='float')
# f3 = np.loadtxt(sys.argv[3], dtype='float')
# f4 = np.loadtxt(sys.argv[4], dtype='float')
# f5 = np.loadtxt(sys.argv[5], dtype='float')
# f6 = np.loadtxt(sys.argv[6], dtype='float')
print("data to be averaged:{}".format(len(sys.argv)-1))
if len(sys.argv) == 11:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f5 = np.loadtxt(sys.argv[5], dtype='float')
    f6 = np.loadtxt(sys.argv[6], dtype='float')
    f7 = np.loadtxt(sys.argv[7], dtype='float')
    f8 = np.loadtxt(sys.argv[8], dtype='float')
    f9 = np.loadtxt(sys.argv[9], dtype='float')
    f10 = np.loadtxt(sys.argv[10], dtype='float')
    f_ave = (f1+f2+f3+f4+f5+f6+f7+f8+f9+f10)/10.0
if len(sys.argv) == 10:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f5 = np.loadtxt(sys.argv[5], dtype='float')
    f6 = np.loadtxt(sys.argv[6], dtype='float')
    f7 = np.loadtxt(sys.argv[7], dtype='float')
    f8 = np.loadtxt(sys.argv[8], dtype='float')
    f9 = np.loadtxt(sys.argv[9], dtype='float')
    f_ave = (f1+f2+f3+f4+f5+f6+f7+f8+f9+f9)/9.0
elif len(sys.argv) == 6:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f5 = np.loadtxt(sys.argv[5], dtype='float')
    f_ave=(f1+f2+f3+f4+f5)/5.0
elif len(sys.argv) == 7:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f5 = np.loadtxt(sys.argv[5], dtype='float')
    f6 = np.loadtxt(sys.argv[6], dtype='float')
    f_ave = (f1+f2+f3+f4+f5+f6)/6.0
elif len(sys.argv) == 8:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f5 = np.loadtxt(sys.argv[5], dtype='float')
    f6 = np.loadtxt(sys.argv[6], dtype='float')
    f7 = np.loadtxt(sys.argv[6], dtype='float')
    f_ave = (f1+f2+f3+f4+f5+f6+f7)/7.0
elif len(sys.argv) == 9:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f5 = np.loadtxt(sys.argv[5], dtype='float')
    f6 = np.loadtxt(sys.argv[6], dtype='float')
    f7 = np.loadtxt(sys.argv[6], dtype='float')
    f8 = np.loadtxt(sys.argv[6], dtype='float')
    f_ave = (f1+f2+f3+f4+f5+f6+f7+f8)/8.0
elif len(sys.argv) == 5:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f4 = np.loadtxt(sys.argv[4], dtype='float')
    f_ave = (f1+f2+f3+f4)/4.0
elif len(sys.argv) == 4:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f3 = np.loadtxt(sys.argv[3], dtype='float')
    f_ave = (f1+f2+f3)/3.0
# elif len(sys.argv) == 3:
#     f1 = np.loadtxt(sys.argv[1], dtype='float')
#     f2 = np.loadtxt(sys.argv[2], dtype='float')
#     f3 = np.loadtxt(sys.argv[3], dtype='float')
#     f_ave = (f1 + f2 + f3)/3.0
else:
    f1 = np.loadtxt(sys.argv[1], dtype='float')
    f2 = np.loadtxt(sys.argv[2], dtype='float')
    f_ave = (f1 + f2)/2.0
# print(f_ave) 
# f_ave.tofile("average.dat")
np.savetxt("average.dat", f_ave)
