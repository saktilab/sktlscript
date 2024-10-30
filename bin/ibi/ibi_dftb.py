#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import math
from math import factorial
from scipy.interpolate import CubicSpline
from rep_spline import RepulsivePotenial
import os
import sys
import pathlib
import argparse

parser = argparse.ArgumentParser(description='IBI')

parser.add_argument('-p', '--prefix', type=str)
parser.add_argument('-i', '--input_skfile', type=str, required=True)
parser.add_argument('-o', '--output_skfile', type=str)
parser.add_argument('-r', '--reference_rdf', type=str, required=True)
parser.add_argument('-c', '--calculated_rdf', type=str, required=True)
parser.add_argument('-f', '--force', action='store_true')
parser.add_argument('-t', '--temperature', type=float, default=300.0)
parser.add_argument('-s', '--start', type=float)
parser.add_argument('-e', '--end', type=float)
opt = parser.parse_args(sys.argv[1:])

temperature_in_K = opt.temperature
ibi_cutoff_in_angstrom = 5.0
ibi_starting_in_angstrom = 2.3
if opt.start :
    ibi_starting_in_angstrom = opt.start
if opt.end:
    ibi_cutoff_in_angstrom = opt.end

output_prefix = 'outputs'
if opt.prefix :
    if pathlib.Path(opt.prefix).exists():
        if opt.force == False:
            raise RuntimeError('Prefix exists')
    output_prefix = opt.prefix


reference_rdf_filename = opt.reference_rdf  #'../cp2k/rdf_oo_py.dat'
current_rdf_filename = opt.calculated_rdf   #'dftb/iter_0_big/rdf_oo_part.dat'
current_sk_filename = opt.input_skfile      #'dftb/iter_0_big/O-O.skf'

os.makedirs(output_prefix, exist_ok=True)
root_path = pathlib.Path(output_prefix)

if not root_path.is_dir():
    raise RuntimeError(f'{output_prefix} is not a folder.')

output_skfile = 'ibi_output.skf'
if opt.output_skfile:
    output_skfile = opt.output_skfile
next_sk_filename = root_path.joinpath(output_skfile)



#figurenames
fig_rdf_name = root_path.joinpath('plot_rdf.svg')
fig_smoothed_rdf_name = root_path.joinpath('plot_smoothed_rdf.svg')
fig_spline_rdf_name = root_path.joinpath('plot_spline_rdf.svg')
fig_deltaf_name = root_path.joinpath('plot_delta_f.svg')
fig_unsmooth_rep_name = root_path.joinpath('plot_unsmooth_rep.svg')
fig_final_rep_name = root_path.joinpath('plot_final_rep.svg')

#constants
bohr_to_angstrom = 0.529177
boltzmann_in_hartree = 3.166811429e-6

#options
show_figures = False
show_intermediate_figures = False


def loadrdf(filename):
    with open(filename, 'r') as f:
        x, y = np.loadtxt(f, usecols=(0, 1), unpack=True)
    return (x, y)
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    if not (isinstance(window_size, int) and isinstance(order, int)):
        raise ValueError("window_size and order must be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError('window_size must be a positive odd number')
    if window_size < order + 2:
        raise TypeError('window_size is too small for the polynomials order')

    order_range = range(order+1)
    half_window = (window_size - 1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window,
                                                           half_window +1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

def write_final_sk(previous_skfile, filename, ori_skfile, current_splines, real_cutoff, resolution=0.01):
    #Fit the first derivatives
    
    nknots = int((real_cutoff - ori_skfile.knots[0])/resolution)+1
    final_knots = np.linspace(ori_skfile.knots[0], real_cutoff, nknots)
    pot_values = [current_splines(r, 0) for r in final_knots]
    first_derivs = [current_splines(r, 1) for r in final_knots]
    second_deriv = ori_skfile.eval(final_knots[0], 2)
    first_derivs_splines = CubicSpline(final_knots, first_derivs, bc_type=((1, second_deriv), (1, 0.0)))

    a4, a3, a2, a1 = first_derivs_splines.c

    r = ori_skfile.knots[0]
    new_expA = -first_derivs_splines(r,1)/first_derivs_splines(r) 
    new_expB = math.log(-first_derivs_splines(r)/new_expA) + new_expA*r
    new_expC = current_splines(r, 0) - math.exp(-new_expA*r+new_expB)


    with open(filename, 'w') as f, open(previous_skfile, 'r') as fin:
        for line in fin:
            if 'Spline' in line:
                break
            print(line, end='', file=f)
            
        print('Spline4', file=f)
        print('{} {:<7.4f}'.format(len(a1), real_cutoff), file=f)
        print(new_expA, new_expB, new_expC, file=f)
    
        for i in range(len(a1)):
            print('{:<7.4f} {:<7.4f} {:18.10E} {:18.10E} {:18.10E} {:18.10E} {:18.10E}'.format(
            first_derivs_splines.x[i], first_derivs_splines.x[i+1], pot_values[i], a1[i], a2[i]/2.0, a3[i]/3.0, a4[i]/4.0
            ), file=f)


target_rdf = loadrdf(reference_rdf_filename)
calc_rdf = loadrdf(current_rdf_filename)

plt.plot(target_rdf[0], target_rdf[1], 'b-', label='Ref. RDF.')
plt.plot(calc_rdf[0], calc_rdf[1], 'r-', label='DFTB RDF.')
plt.title('RDF')
plt.legend()

if (show_figures and show_intermediate_figures):
    plt.show()
    plt.savefig(fig_rdf_name, format='svg')
    plt.clf()
else:
    plt.savefig(fig_rdf_name, format='svg')
    plt.clf()

target_ys = savitzky_golay(target_rdf[1], 19, 3)
plt.plot(target_rdf[0], target_rdf[1], 'b-', label='Ref. RDF.')
plt.plot(target_rdf[0], target_ys, 'r-', label='Smoothed Ref. RDF.')
calc_ys = savitzky_golay(calc_rdf[1], 19, 3)
plt.plot(calc_rdf[0], calc_rdf[1], 'b--', label='DFTB RDF.')
plt.plot(calc_rdf[0], calc_ys, 'r--', label='Smoothed DFTB RDF.')
plt.title('Smoothed RDF')
plt.legend()

if (show_figures and show_intermediate_figures):
    plt.show()
    plt.savefig(fig_smoothed_rdf_name, format='svg')
    plt.clf()
else:
    plt.savefig(fig_smoothed_rdf_name, format='svg')
    plt.clf()


spline_target = CubicSpline(target_rdf[0], target_ys)
spline_calc = CubicSpline(calc_rdf[0], calc_ys)
plt.plot(calc_rdf[0], calc_rdf[1], 'b-', label='Ref. RDF.')
plt.plot(calc_rdf[0], spline_calc(calc_rdf[0]), 'r--', label='Spline Interpolated Ref. RDF.')
plt.title('Spline-interpolated Smoothed RDF')
plt.legend()

if (show_figures and show_intermediate_figures):
    plt.show()
    plt.savefig(fig_spline_rdf_name, format='svg')
    plt.clf()
else:
    plt.savefig(fig_spline_rdf_name, format='svg')
    plt.clf()


#section of IBI

kbt = boltzmann_in_hartree*temperature_in_K

cutoff = ibi_cutoff_in_angstrom/bohr_to_angstrom
ibi_starting = ibi_starting_in_angstrom/bohr_to_angstrom

ibi_grids = np.linspace(ibi_starting, cutoff, 200)

def ibi(r):  # r in bohrm splines are in angstrom
    try:
        g_ref = spline_target(r*bohr_to_angstrom)
        g_calc = spline_calc(r*bohr_to_angstrom)
        delta_f = kbt*math.log(g_calc/g_ref)
        return delta_f
    except:
        return 0.0

fs = [ibi(r) for r in ibi_grids]
    
plt.plot(ibi_grids, fs, 'r-')
plt.title('Delta F')

if (show_figures and show_intermediate_figures):
    plt.show()
    plt.savefig(fig_deltaf_name, format='svg')
    plt.clf()
else:
    plt.savefig(fig_deltaf_name, format='svg')
    plt.clf()


skfile = RepulsivePotenial.from_file(current_sk_filename)


#smooth the final rep
smooth_windows = 2.0
smoothing_grids = np.linspace(ibi_starting-smooth_windows, ibi_starting+smooth_windows, 41).tolist()
final_rep_wo_smooth = []
ori_rep = []
fs = []
for r in smoothing_grids:
    if r < ibi_starting:
        fs.append(0.0)
        final_rep_wo_smooth.append(skfile.eval(r))
    else:
        fs.append(ibi(r))
        final_rep_wo_smooth.append(skfile.eval(r)+ibi(r))
    ori_rep.append(skfile.eval(r))

final_rep_smooth = savitzky_golay(np.array(final_rep_wo_smooth), 9, 5)
final_rep_smooth = savitzky_golay(np.array(final_rep_smooth), 5, 3)

plt.plot(smoothing_grids, ori_rep, 'r-', label='Original Rep.')
plt.plot(smoothing_grids, final_rep_wo_smooth, 'bo-', markersize=3, label='Next Rep.')
plt.plot(smoothing_grids, final_rep_smooth, 'go-', markersize=2, label='Smoothed Next Rep.')


plt.title('Repulsive Potential')
plt.xlim(ibi_starting-smooth_windows-1.0, smoothing_grids[-1]+1.0)
plt.ylim(-0.002, 0.005)
plt.legend()

if (show_figures):
    plt.savefig(fig_unsmooth_rep_name, format='svg')
    plt.show()
    plt.clf()
else:
    plt.savefig(fig_unsmooth_rep_name, format='svg')
    plt.clf()


#fitting the final rep
fitting_windows = 1.0
fitting_grids = np.linspace(skfile.knots[0], cutoff, 200).tolist()

final_grids = []
final_pot = []
for r in fitting_grids:
    if r < ibi_starting-fitting_windows:
        final_grids.append(r)
        final_pot.append(skfile.eval(r))
    else:
        break

for r, v in zip(smoothing_grids, final_rep_smooth):
    if r >= ibi_starting-fitting_windows and r < ibi_starting+fitting_windows:
        final_grids.append(r)
        final_pot.append(v)
for r in fitting_grids:
    if (r>=ibi_starting+fitting_windows):
        final_grids.append(r)
        final_pot.append(skfile.eval(r)+ibi(r))

real_cutoff = cutoff+0.5
final_grids.append(real_cutoff)
final_pot.append(0.0)


first_deriv = skfile.eval(final_grids[0], 1)
final_splines = CubicSpline(final_grids, final_pot, bc_type=((1, first_deriv), (1, 0.0)))

ploting_grids = np.linspace(skfile.knots[0], real_cutoff, 1000)

plt.plot(smoothing_grids, ori_rep, 'r-', label='Original Rep.')
plt.plot(final_grids, final_pot, 'bo-', markersize=3, label='Next Rep.')
plt.plot(ploting_grids, final_splines(ploting_grids), 'go-', markersize=2, label='Smoothed Next Rep.')

plt.title('Repulsive Potential')
plt.xlim(ibi_starting-1.0, real_cutoff)
plt.ylim(-0.002, 0.005)
plt.legend()

if (show_figures):
    plt.savefig(fig_final_rep_name, format='svg')
    plt.show()
else:
    plt.savefig(fig_final_rep_name, format='svg')
    plt.clf()


write_final_sk(current_sk_filename, next_sk_filename, skfile, final_splines, real_cutoff)
