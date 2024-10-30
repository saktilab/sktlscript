#!/usr/bin/env python3
import scipy.stats  as stats
import sys
import numpy as np
import matplotlib.pyplot as plt  
from scipy.stats import gaussian_kde
from matplotlib import rc

fout = open("delta.dat", "w")
R = 8.314
T = float(sys.argv[2])
with open(sys.argv[1], "r") as f:
    f.readline()
    deltas = []
    force = []
    delta_coord = []
    for line in f:
        arr = line.split()
        delta = float(arr[3])
        deltas.append(delta)    
    
    # mindelta = np.min(deltas)
    # maxdelta = np.max(deltas)
    # mean = np.mean(deltas)
    # std = np.std(deltas)
    # prob = np.histogram(deltas, bins=100, range=[mindelta, maxdelta], density=True)
    # bins, probs = prob
    density = gaussian_kde(deltas)
    xs = np.linspace(0,max(deltas),200)
    density._norm_factor = 1
    density.covariance_factor= lambda : .1
    density._compute_covariance()
    P = density(xs)/sum(density(xs))

    deltaF = -R*T*np.log(P)/1000/4.184 
    deltaF = deltaF - min(deltaF)
    arr = np.column_stack([xs, deltaF])
    # print(('{:15.4f} ' ' {:15.4f}').format(xs, deltaF),file=fout)
    np.savetxt(fout, arr)
    # plt.plot(xs, deltaF)
    # plt.show() )
# plt.figure(figsize=(50, 0))

###########START TO PLOT#####
plt.figure(dpi=300)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

ax = plt.subplot(111)  
ax.spines["top"].set_visible(True)  
ax.spines["right"].set_visible(True) 

# ax.get_xaxis().tick_bottom()  
# ax.get_yaxis().tick_left() 
plt.tick_params(direction='in')
plt.ylim(0, 5)
plt.xlim(0,0.7)
# plt.xticks(range(0.0, 0.9, 0.1), fontsize=14)  
plt.yticks(np.arange(0.5,5, step=0.5), fontsize=14)
plt.xticks(fontsize=14)
plt.ylabel(r'$\Delta F$ [kcal/mol]', fontsize=18)
plt.xlabel(r'$\delta$ [\AA]', fontsize=18)
plt.plot(xs, deltaF, color="blue", lw=2)

plt.savefig("pmf.pdf", bbox_inches="tight");  

    
    