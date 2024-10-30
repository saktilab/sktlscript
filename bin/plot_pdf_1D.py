#!/usr/bin/env python3
import scipy.stats  as stats
import sys
import numpy as np
import matplotlib.pyplot as plt  
from scipy.stats import gaussian_kde
from matplotlib import rc

R = 8.314
with open(sys.argv[1], "r") as f:
    f.readline()
    charges = []
    for line in f:
        arr = line.split()
        charge = float(arr[1])
        charges.append(charge)    
    
    # mindelta = np.min(deltas)
    # maxdelta = np.max(deltas)
    # mean = np.mean(deltas)
    # std = np.std(deltas)
    # prob = np.histogram(deltas, bins=100, range=[mindelta, maxdelta], density=True)
    # bins, probs = prob
    density = gaussian_kde(charges)
    xs = np.linspace(0,max(charges),sys.argv[2])
    density._norm_factor = 1
    density.covariance_factor= lambda : .1
    density._compute_covariance()
    P = density(xs)/sum(density(xs))
    data = np.column_stack([xs, P])
    np.savetxt("pdf.dat", data)
###########START TO PLOT#####
plt.figure(dpi=300)
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

ax = plt.subplot(111)  
ax.spines["top"].set_visible(True)  
ax.spines["right"].set_visible(True) 

# ax.get_xaxis().tick_bottom()  
# ax.get_yaxis().tick_left() 
plt.tick_params(direction='in')
# plt.ylim(0, 5)
# plt.xlim(0,0.7)
# plt.xticks(range(0.0, 0.9, 0.1), fontsize=14)  
plt.yticks(np.arange(0.5,5, step=0.5), fontsize=14)
plt.xticks(fontsize=14)
# plt.ylabel('Probability density', fontsize=18)
# plt.xlabel('Charge of Ce', fontsize=18)
plt.plot(xs, P, color="blue", lw=2)

plt.savefig("pdf.pdf", bbox_inches="tight");  

    
    