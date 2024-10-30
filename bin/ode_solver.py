#!/usr/bin/env python3
import numpy as np
from scipy.integrate import odeint

def rhs(y,t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14):
	r1 = k1*y[0]
	r2 = k2*y[1]
	r3 = k3*y[1]
	r4 = k4*y[2]
	r5 = k5*y[2]
	r6 = k6*y[3]
	r7 = k7*y[3]
	r8 = k8*y[4]*y[5]
	r9 = k9*y[4]*y[6]
	r10 = k10*y[7]
	r11 = k11*y[7]
	r12 = k12*y[8]
	r13 = k13*y[8]
	r14 = k14*y[9]
	# to return the rate equations of EtOH, B, C, D, E, F, Ace, FAL, G, H, I
	return [r2-r1,r1-r2-r3+r4,r3-r4-r5+r6,r5-r6-r7+r8,r7-r8-r9+r10,r7-r8-r9,r10-r9,r9-r10-r11+r12,r11-r12-r13+r14,r13-r14]

time = np.linspace(0,100)
k_vals = 8.13E+11, 9.72E+12, 3.37E+00, 2.30E+11, 1.63E+10, 1.52E+04, 3.82E+15, 8.08E+08, 6.06E+09, 4.11E+14, 3.46E-03, 1.24E-01, 6.93E+11, 2.51E+08
y0 = [0.2,0,0,0,0,0,0,0,0,0]
yout = odeint(rhs,y0,time,k_vals)
import matplotlib.pyplot as plt


plt.plot(time,yout)
with open("timecourse.dat",'w') as f:
	print(time,yout[0:10],file=f)
_ = plt.legend(['EtOH','C','D','E','F','Ace','FAL','G','H','I'])
# plt.ylim(0,)
plt.show()
