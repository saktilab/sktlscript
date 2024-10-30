#!/usr/bin/env python3

import sys
import random
import numpy as np

class PolynomialRepulsive:
    def __init__(self, cutoff, coeffs, minimal_order=3):
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
    

