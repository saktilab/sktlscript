#!/usr/bin/env python3
from pymatgen.io.vasp import Vasprun
# from pymatgen.electronic_structure.plotter import DosPlotter
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import matplotlib.pylab  as pyl
import sys
# import sys
energy_dftb = []
density_dftb = []
# with open("dos_total.dat", "r") as f:
#     for line in f:
#         arr = line.split()
#         energy_dftb.append(arr[0])
#         density_dftb.append(arr[1])

fontlabel=18
fonttick=14

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-4, 6)}

v = Vasprun('vasprun.xml')
tdos = v.tdos

# densities = tdos.get_densities()
densities = tdos.get_smeared_densities(float(sys.argv[1]))
# densities = densities.values
energies = tdos.energies
shifted = []
for energy in energies:
    shifted.append(energy-tdos.efermi)
# ala = []
smeared_densities = []
densities = list(densities.values())
for i in densities:
    for x in i:
        smeared_densities.append(x)
for energy,density in zip(shifted,smeared_densities):
    print(energy, density)
# com_dens = v.complete_dos.get_smeared_densities(float(sys.argv[1]))
# com_energies = v.complete_dos.energies
# shifted_com = []
# for energy in com_energies:
#     shifted_com.append(energy - v.complete_dos.efermi)
# Ce_dos = v.complete_dos.get_element_spd_dos("Ce")

# print(shifted,densities)
# xlabel = "Energy [eV]"
# ylabel = "DOS [arb. units]"

# plt.rcParams.update(params)
# plt.subplots_adjust(left=0.21, right=0.93,
#                     bottom=0.16, top=0.88,
#                     wspace=None, hspace=None)

# yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)                           
# figure = plt.figure(1, figsize=(8,5.5))
# ax     = figure.add_subplot(111)
# ax.text(0.5,-0.15,xlabel,size=fontlabel,
#         transform=ax.transAxes,ha='center',va='center')
# ax.text(-0.15,0.5,ylabel,size=fontlabel,
#         transform=ax.transAxes,ha='center',va='center',rotation=90)
# for line in ax.get_xticklines() + ax.get_yticklines():
#     line.set_markersize(6)
#     line.set_markeredgewidth(2)
# plt.xticks(size=fonttick)
# plt.yticks(size=fonttick)
# # pyl.grid(True)
# # plt.ylim(0,110)
# # plt.xlim(-5,3)
# # plt.ylim(0,40)
# # plt.plot(shifted,densities,'b-',label="PBE+U")
# # plt.plot(energy_dftb, density_dftb, 'b-', label="DFTB")
# # plt.legend(loc=9,borderaxespad=.8,numpoints=1)

# plt.savefig('tdos.png',orientation='portrait',format='png',dpi=300,bbox_inches = 'tight',
#     pad_inches = 0.01)
