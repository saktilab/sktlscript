#!/usr/bin/env python3


import sys
import matplotlib
import argparse
import glob
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import matplotlib
#rcParams['text.usetex'] = True
#rcParams['text.latex.unicode'] = True
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['axes.linewidth'] = 2
font = {'family' : 'serif',
        'serif': 'Times',
        'size': 16}

#matplotlib.rc('font', **font)

kpaths = {
    'bcc': [r'$\Gamma$', 'H', 'N', r'$\Gamma$', 'P', 'H,P', 'N'],
    'fcc': [r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L', 'U', 'W', 'L', 'K,U', 'X'],
    'hcp': [r'$\Gamma$', 'M', 'K', r'$\Gamma$', 'A', 'L', 'H', 'A,L', 'M,K', 'H'],
    'sc' : [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R', 'X,M', 'R']
}

my_colors = ['k','r','b','g', 'c', 'm'] 

parser = argparse.ArgumentParser(description='plot bandstructure')
# #parser.add_argument('-n', default=0, type=int)
#parser.add_argument('-c', metavar='cell_type', type=str)
parser.add_argument('-yx', '--ymax', type=float)
parser.add_argument('-yi', '--ymin', type=float)
parser.add_argument('-c', '--cell_type', type=str)
parser.add_argument('-g', '--geometry', type=str)
parser.add_argument('-t', '--title', type=str)
parser.add_argument('plots', metavar='plots', nargs='+', type=str)
opt = parser.parse_args()

ymax = 10.0
ymin = -12.0
if opt.ymax :
    ymax = opt.ymax
if opt.ymin :
    ymin = opt.ymin

print (opt)

def loadnxy(filenames):
    filename = glob.glob(filenames)[0]
    temp_lines = []
    with open(filename, 'r') as f:
        for line in f:
            arr = line.split()
            line_data = list(map(float, arr[1:]))
            temp_lines.append(line_data)
    return np.transpose(temp_lines)         

#cell_type = sys.argv[1]
#plots = sys.argv[2:]
cell_type = opt.cell_type

plots = opt.plots

plt.axhline(0, color='black', linestyle='--')

ticks_pos = []

legend_plots = []
for ind, plot in enumerate(plots):
    filename, my_label,  = plot.split(':')

    
    data = loadnxy(filename)
    nkpt = np.size(data,1)

    xs = list(range(nkpt))
    for pi, band in enumerate(data):
        handle, = plt.plot(xs, band, color=my_colors[ind], label=my_label)
        if (pi==0):
            legend_plots.append(handle)

    if (cell_type is not None and len(ticks_pos)==0):
        num_label = len(kpaths[cell_type])
        print (nkpt , num_label-1)
        if (cell_type in kpaths and (nkpt % (num_label-1) == 0)):
            interval = nkpt / (num_label-1)
            for i in range(len(kpaths[cell_type])):
                ticks_pos.append(i*interval)  

if (cell_type is not None):
    print (ticks_pos, kpaths[cell_type])
    plt.xticks(ticks_pos, kpaths[cell_type])

plt.ylabel('Energy [eV]')
plt.ylim(ymax=ymax, ymin=ymin) 
plt.legend(handles=legend_plots)

if (opt.title is not None):
    plt.title(opt.title,y=1.02)

ax = plt.gca()
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.xaxis.grid(True, which='major')
ax.xaxis.set_tick_params(width=2, length=5)
ax.yaxis.set_tick_params(width=2, length=5)
# ax.grid(True)
plt.savefig('bs.eps', format='eps', dpi=600)
plt.show()





