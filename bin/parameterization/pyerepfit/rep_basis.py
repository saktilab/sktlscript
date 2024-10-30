#!/usr/bin/env python3
import numpy as np
from rep_spline import RepulsivePotenial
from scipy.interpolate import CubicSpline

class BasisRepulsiveModel:
    def __init__(self, cutoff, minimal_order = 3, order = 9):
        self.cutoff = cutoff
        self.minimal_order = minimal_order
        self.order = order
        self.num_variables = self.order - self.minimal_order + 1
    def get_row(self, r, derv=0):
        res = np.zeros(self.num_variables)
        for i in range(self.num_variables):
            if (derv == 0):
                res[i] =  (self.cutoff-r)**(i+self.minimal_order)
            elif (derv == 1):
                res[i] =  -1.0*(i+self.minimal_order)*(self.cutoff-r)**(i+self.minimal_order-1)
            elif (derv == 2):
                res[i] =  (i+self.minimal_order)*(i+self.minimal_order-1)*(self.cutoff-r)**(i+self.minimal_order-2)
            else:
                raise ValueError('Not valid derivative')
        return res
    def description(self):
        return "Polynomial-type repulsive potential (R-R_cutoff)^n"
    
    def __str__(self):
        return "Minimal and maximal orders: {}, {}.\nCutoff: {}".format(self.minimal_order, self.order, self.cutoff)

class PolynomialRepulsive:
    def __init__(self, name, cutoff, coeffs, minimal_order=3):
        self.name = name
        self.cutoff = cutoff
        self.coeffs = coeffs
        self.minimal_order = minimal_order
    
    def eval(self, distance: float, derivs=0):
        if (distance >= self.cutoff):
            return 0.0
        func = 0.0
        if (derivs == 0):
            for i, c in enumerate(self.coeffs):
                func += c*(self.cutoff-distance)**(self.minimal_order+i) 
        elif (derivs==1):
            for i, c in enumerate(self.coeffs):
                func += -1.0*(self.minimal_order+i)*c*(self.cutoff-distance)**(self.minimal_order+i-1) 
        elif (derivs==2):
            for i, c in enumerate(self.coeffs):
                func += (self.minimal_order+i)*(self.minimal_order+i-1)*c*(self.cutoff-distance)**(self.minimal_order+i-2) 
        else:
            raise ValueError('Derivative out of range')
        return func

    def to_spline4(self, interval=0.05, start=1.0):
        num = (self.cutoff - start)/interval + 1
        xs = np.linspace(1.0, self.cutoff, num)
        ys = [self.eval(r, 1) for r in xs]
        d1 = self.eval(xs[0], 2)
        y0s = [self.eval(r, 0) for r in xs]
        # print(self.name)
        # with open('{}.rep'.format(self.name), 'w') as f:
        #     for x, y in zip(xs, y0s):
        #         print('{:10.5f} {:20.12f}'.format(x, y), file=f)
        ala = CubicSpline(xs, ys, bc_type=((1, d1), (1, 0.0)))
        spline_coeffs = []
        for i in range(len(ala.c[0])):
            l_coeffs = [ self.eval(xs[i], 0), ala.c[3][i], ala.c[2][i]/2.0, ala.c[1][i]/3.0, ala.c[0][i]/4.0]
            spline_coeffs.append(l_coeffs)
        spl_rep = RepulsivePotenial(self.name, xs, spline_coeffs)
        return spl_rep
