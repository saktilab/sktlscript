#!/usr/bin/env python3
from pymatgen.io.vasp import Vasprun
# from pymatgen.electronic_structure.plotter import DosPlotter
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import matplotlib.pylab  as pyl
import sys
from pymatgen.electronic_structure.core import Spin, Orbital
from numpy import array as npa
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
cdos = v.complete_dos

# densities = tdos.get_densities()
densities = cdos.get_smeared_densities(float(sys.argv[1]))
# densities = densities.values
cbm_vbm = [x for x in cdos.get_cbm_vbm()]
band_gap = cbm_vbm[0]-cbm_vbm[1]
print("Direct band gap = {:.2f} eV".format(band_gap))

energies = cdos.energies
shifted = []
for energy in energies:
    shifted.append(energy-cbm_vbm[0])
#Collect densities 
density_Ce_s = list(v.pdos[0][Orbital.s][Spin.up])

density_Ce_p = list(npa(v.pdos[0][Orbital.px][Spin.up])+npa(v.pdos[0][Orbital.py][Spin.up])+npa(v.pdos[0][Orbital.pz][Spin.up]))

density_Ce_f = list(npa(v.pdos[0][Orbital.f0][Spin.up])+npa(v.pdos[0][Orbital.f1][Spin.up])+npa(v.pdos[0][Orbital.f2][Spin.up])+npa(v.pdos[0][Orbital.f3][Spin.up])+npa(v.pdos[0][Orbital.f_1][Spin.up])+npa(v.pdos[0][Orbital.f_2][Spin.up])+npa(v.pdos[0][Orbital.f_3][Spin.up]))

density_Ce_d = list(npa(v.pdos[0][Orbital.dx2][Spin.up])+npa(v.pdos[0][Orbital.dxy][Spin.up])+npa(v.pdos[0][Orbital.dxz][Spin.up])+npa(v.pdos[0][Orbital.dyz][Spin.up])+npa(v.pdos[0][Orbital.dz2][Spin.up]))

density_O_s = list(v.pdos[1][Orbital.s][Spin.up])

density_O_p = list(npa(v.pdos[1][Orbital.px][Spin.up])+npa(v.pdos[1][Orbital.py][Spin.up])+npa(v.pdos[1][Orbital.pz][Spin.up]))

if (sys.argv[2]=="spin"):
    density_Ce_s_up = list(v.pdos[0][Orbital.s][Spin.up])
    density_Ce_s_down = list(v.pdos[0][Orbital.s][Spin.down])

    density_Ce_p_up = list(npa(v.pdos[0][Orbital.px][Spin.up])+npa(v.pdos[0][Orbital.py][Spin.up])+npa(v.pdos[0][Orbital.pz][Spin.up]))
    density_Ce_p_down = list(npa(v.pdos[0][Orbital.px][Spin.up])+npa(v.pdos[0][Orbital.py][Spin.up])+npa(v.pdos[0][Orbital.pz][Spin.down]))

    density_Ce_f_up = list(npa(v.pdos[0][Orbital.f0][Spin.up])+npa(v.pdos[0][Orbital.f1][Spin.up])+npa(v.pdos[0][Orbital.f2][Spin.up])+npa(v.pdos[0][Orbital.f3][Spin.up])+npa(v.pdos[0][Orbital.f_1][Spin.up])+npa(v.pdos[0][Orbital.f_2][Spin.up])+npa(v.pdos[0][Orbital.f_3][Spin.up]))
    
    density_Ce_f_down = list(npa(v.pdos[0][Orbital.f0][Spin.up])+npa(v.pdos[0][Orbital.f1][Spin.up])+npa(v.pdos[0][Orbital.f2][Spin.up])+npa(v.pdos[0][Orbital.f3][Spin.up])+npa(v.pdos[0][Orbital.f_1][Spin.up])+npa(v.pdos[0][Orbital.f_2][Spin.up])+npa(v.pdos[0][Orbital.f_3][Spin.down]))

    density_Ce_d_up = list(npa(v.pdos[0][Orbital.dx2][Spin.up])+npa(v.pdos[0][Orbital.dxy][Spin.up])+npa(v.pdos[0][Orbital.dxz][Spin.up])+npa(v.pdos[0][Orbital.dyz][Spin.up])+npa(v.pdos[0][Orbital.dz2][Spin.up]))
    
    density_Ce_d_down = list(npa(v.pdos[0][Orbital.dx2][Spin.up])+npa(v.pdos[0][Orbital.dxy][Spin.up])+npa(v.pdos[0][Orbital.dxz][Spin.up])+npa(v.pdos[0][Orbital.dyz][Spin.up])+npa(v.pdos[0][Orbital.dz2][Spin.down]))
    

    density_O_s_up = list(v.pdos[1][Orbital.s][Spin.up])
    density_O_s_down = list(v.pdos[1][Orbital.s][Spin.down])
    
    density_O_p_up = list(npa(v.pdos[1][Orbital.px][Spin.up])+npa(v.pdos[1][Orbital.py][Spin.up])+npa(v.pdos[1][Orbital.pz][Spin.up]))
    density_O_p_down = list(npa(v.pdos[1][Orbital.px][Spin.up])+npa(v.pdos[1][Orbital.py][Spin.up])+npa(v.pdos[1][Orbital.pz][Spin.down]))
    
    with open("dos_ce_s_up.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_s_up):
            print(energy, density, file=fout)

    with open("dos_ce_p_up.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_p_up):
            print(energy, density, file=fout)

    with open("dos_ce_f_up.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_f_up):
            print(energy, density, file=fout)

    with open("dos_ce_d_up.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_d_up):
            print(energy, density, file=fout)

    with open("dos_o_s_up.dat", "w") as fout:
        for energy, density in zip(shifted, density_O_s_up):
            print(energy, density, file=fout)

    with open("dos_o_p_up.dat", "w") as fout:
        for energy, density in zip(shifted, density_O_p_up):
            print(energy, density, file=fout)

    with open("dos_ce_s_down.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_s_down):
            print(energy, -density, file=fout)

    with open("dos_ce_p_down.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_p_down):
            print(energy, -density, file=fout)

    with open("dos_ce_f_down.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_f_down):
            print(energy, -density, file=fout)

    with open("dos_ce_d_down.dat", "w") as fout:
        for energy, density in zip(shifted, density_Ce_d_down):
            print(energy, -density, file=fout)

    with open("dos_o_s_down.dat", "w") as fout:
        for energy, density in zip(shifted, density_O_s_down):
            print(energy, -density, file=fout)

    with open("dos_o_p_down.dat", "w") as fout:
        for energy, density in zip(shifted, density_O_p_down):
            print(energy, -density, file=fout)
########################################
with open("dos_ce_s.dat", "w") as fout:
    for energy, density in zip(shifted, density_Ce_s):
        print(energy, density, file=fout)

with open("dos_ce_p.dat", "w") as fout:
    for energy, density in zip(shifted, density_Ce_p):
        print(energy, density, file=fout)

with open("dos_ce_f.dat", "w") as fout:
    for energy, density in zip(shifted, density_Ce_f):
        print(energy, density, file=fout)

with open("dos_ce_d.dat", "w") as fout:
    for energy, density in zip(shifted, density_Ce_d):
        print(energy, density, file=fout)

with open("dos_o_s.dat", "w") as fout:
    for energy, density in zip(shifted, density_O_s):
        print(energy, density, file=fout)

with open("dos_o_p.dat", "w") as fout:
    for energy, density in zip(shifted, density_O_p):
        print(energy, density, file=fout)



# smeared_densities = []
# densities = list(densities.values())



# for i in densities:
#     for x in i:
#         smeared_densities.append(x)
# for energy,density in zip(shifted,smeared_densities):
#     print(energy, density)
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
