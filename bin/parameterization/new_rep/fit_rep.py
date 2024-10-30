#!/usr/bin/env python3

import sys
import random
import numpy as np
import spline_rep
import math
import pygmo as pg
import new_rep
import os

import scipy.optimize

import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline
class fitness_function:

    def __init__(self, spline_rep, cutoff, max_order=9, minimal_order=3, start=1.0):
        self.max_order = max_order
        self.spline_rep = spline_rep
        self.cutoff = cutoff
        self.minimal_order = minimal_order
        self.start = start

    def fitness(self, x):
        rep = new_rep.PolynomialRepulsive(self.cutoff, x, self.minimal_order)
        error = 0.0
        # xs = np.linspace(self.spline_rep.knots[0], max(self.cutoff, self.spline_rep.knots[-1]), num=400)
        xs = np.linspace(self.start, max(self.cutoff, self.spline_rep.knots[-1]), num=400)
        for r in xs:
            v = rep.eval(r)
            v_ref = self.spline_rep.eval(r)
            error += (v-v_ref)*(v-v_ref)
        error /= float(len(xs)) 
        final_error = math.sqrt(error)
        return [final_error]

    def get_bounds(self):
        num = self.max_order-3+1
        return ([-1]*num,[1]*num)


input_rep = sys.argv[1]
minimal_order = 3
maximal_order = 12

spl_rep = spline_rep.RepulsivePotenial.from_file(input_rep)
f_func = fitness_function(spl_rep, spl_rep.cutoff, maximal_order, minimal_order)

# algo = pg.core.algorithm(pg.sade(gen=2000, ftol=1.0e-7))
algo = pg.core.algorithm(pg.xnes(gen=2000, ftol=1.0e-7))
algo.set_verbosity(10)
prob = pg.core.problem(f_func)
pop = pg.core.population(prob,50)
pop = algo.evolve(pop)

print('#', pop.champion_x)
print('#', pop.champion_f) 

output = '{}.rep'.format(os.path.basename(input_rep))
name = os.path.basename(input_rep).split('.')[0]
output_fig = '{}.eps'.format(name)

def to_spline4(self, interval=0.05, start=1.0):
    num = (self.cutoff - start)/interval + 1
    xs = np.linspace(1.0, self.cutoff, num)
    ys = [self.eval(r, 1) for r in xs]
    d1 = self.eval(xs[0], 2)
    y0s = [self.eval(r, 0) for r in xs]
    ala = CubicSpline(xs, ys, bc_type=((1, d1), (1, 0.0)))
    spline_coeffs = []
    for i in range(len(ala.c[0])):
        l_coeffs = [ self.eval(xs[i], 0), ala.c[3][i], ala.c[2][i]/2.0, ala.c[1][i]/3.0, ala.c[0][i]/4.0]
        spline_coeffs.append(l_coeffs)
    spl_rep = spline_rep.RepulsivePotenial("", xs, spline_coeffs)
    return spl_rep



plt.title(name)
plt.ylim(-5, 5)
with open(output, 'w') as f:
    rep = new_rep.PolynomialRepulsive(spl_rep.knots[-1], pop.champion_x, minimal_order)
    new_spline_rep = to_spline4(rep, interval=0.05, start=spl_rep.knots[0])
    # print(new_spline_rep)
    styles = ['c--', 'r--', 'g--']
    labels = ['Rep. Potential', '1st derivative', '2nd derivative']
    for i in range(3):
        xs = np.linspace(1.0, spl_rep.knots[-1], num=400)
        ref_ys = []
        ys = []
        for r in xs:
            v = new_spline_rep.eval(r, i)
            v_ref = spl_rep.eval(r, i)
            print(r, v_ref, v, file=f)    
            ys.append(v)
            ref_ys.append(v_ref)
        print('', file=f)
        
        plt.plot(xs, ref_ys, 'k-', linewidth=1)
        plt.plot(xs, ys, styles[i], linewidth=1, label=labels[i])
    # plt.legend([None, 'Rep. Potential', None, '1st derivative', None, '2nd derivative'])
    plt.legend(loc='upper right')
    plt.xlabel('Diatomic Distance [a.u.]')
    plt.ylabel('Energy [a.u.]')
    plt.savefig(output_fig, format='eps')
