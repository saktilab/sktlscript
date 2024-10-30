#!/usr/bin/env python3

import sys
import math
import io
import decimal
import os
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

class Spline5:
    """Fifth-order spline function"""
    
    def __init__(self, r0, r1, a0, a1, a2, a3, a4, a5):
        self.r0 = r0
        self.r1 = r1
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a5 = a5

    def eval(self, r, der=0):
        if (r<self.r0 or r>self.r1):
            raise ValueError('r: {}. Out of range {}-{}.'.format(r, self.r0, self.r1))
        if (der <0 or der>3):
            raise ValueError('der: {}. Out of range 0-3.'.format(der))

        dx = r-self.r0
        if der == 0:
            return self.a0+self.a1*dx+self.a2*dx*dx+self.a3*dx**3+self.a4*dx**4+self.a5*dx**5
        elif der == 1:
            return self.a1+2.0*self.a2*dx+3.0*self.a3*dx*dx+4.0*self.a4*dx**3+5.0*self.a5*dx**4
        elif der == 2:
            return 2.0*self.a2+3.0*2.0*self.a3*dx+4.0*3.0*self.a4*dx*dx+5.0*4.0*self.a5*dx**3
        elif der == 3:
            return 3.0*2.0*self.a3+4.0*3.0*2.0*self.a4*dx+5.0*4.0*3.0*self.a5*dx*dx
    def __str__(self):
        return '{:<12.6f} {:<12.6f} {:20.12E} {:20.12E} {:20.12E} {:20.12E} {:20.12E} {:20.12E}'.format(
            self.r0, self.r1,
            self.a0, self.a1, self.a2, self.a3, self.a4, self.a5
        )

class Spline4:
    """Fourth-order spline function"""
    
    def __init__(self, r0, r1, a0, a1, a2, a3, a4):
        self.r0 = r0
        self.r1 = r1
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4

    def eval(self, r, der=0):
        if (r<self.r0 or r>self.r1):
            raise ValueError('r: {}. Out of range {}-{}.'.format(r, self.r0, self.r1))
        if (der <0 or der>3):
            raise ValueError('der: {}. Out of range 0-3.'.format(der))

        dx = r-self.r0
        if der == 0:
            return self.a0+self.a1*dx+self.a2*dx*dx+self.a3*dx*dx*dx+self.a4*dx*dx*dx*dx
        elif der == 1:
            return self.a1+2.0*self.a2*dx+3.0*self.a3*dx*dx+4.0*self.a4*dx*dx*dx
        elif der == 2:
            return 2.0*self.a2+3.0*2.0*self.a3*dx+4.0*3.0*self.a4*dx*dx
        elif der == 3:
            return 3.0*2.0*self.a3+4.0*3.0*2.0*self.a4*dx
    def __str__(self):
        return '{:<12.6f} {:<12.6f} {:20.12E} {:20.12E} {:20.12E} {:20.12E} {:20.12E}'.format(
            self.r0, self.r1,
            self.a0, self.a1, self.a2, self.a3, self.a4
        )

class Spline3:
    """Cubic spline function"""
    
    def __init__(self, r0, r1, a0, a1, a2, a3):
        self.r0 = r0
        self.r1 = r1
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3

    def eval(self, r, der=0):
        if (r<self.r0 or r>self.r1):
            raise ValueError('r: {}. Out of range {}-{}.'.format(r, self.r0, self.r1))
        if (der <0 or der>2):
            raise ValueError('der: {}. Out of range 0-2.'.format(der))

        dx = r-self.r0
        if der == 0:
            return self.a0+self.a1*dx+self.a2*dx*dx+self.a3*dx*dx*dx
        elif der == 1:
            return self.a1+2.0*self.a2*dx+3.0*self.a3*dx*dx
        elif der == 2:
            return 2.0*self.a2+3.0*2.0*self.a3*dx

    def __str__(self):
        return '{:<12.6f} {:<12.6f} {:20.12E} {:20.12E} {:20.12E} {:20.12E}'.format(
            self.r0, self.r1,
            self.a0, self.a1, self.a2, self.a3
        )


class RepulsivePotenial:
    """Repulsive potential in DFTB"""

    def eval(self, r, der=0):
        if (r < self.splines[0].r0):
            if (der == 0):
                return math.exp(-self.expA*r+self.expB) + self.expC
            elif (der == 1):
                return -self.expA*math.exp(-self.expA*r+self.expB) 
            elif (der == 2):
                return self.expA*self.expA*math.exp(-self.expA*r+self.expB) 
            elif (der == 3):
                return (-self.expA**3)*math.exp(-self.expA*r+self.expB) 
            
        for spline in self.splines:
            if (r >= spline.r0 and r < spline.r1):
                return spline.eval(r, der)
        return 0.0

    def fix_exponential(self):
        found = False
        r = self.knots[0]
        max_r = 2.0
        last_good = None
        while(r < max_r):
            new_expA = -self.eval(r,2)/self.eval(r,1) 
            try:
                new_expB = math.log(-self.eval(r, 1)/new_expA) + new_expA*r
                new_expC = self.eval(r, 0) - math.exp(-new_expA*r+new_expB)
                last_good = (new_expA, new_expB, new_expC)
                if (new_expA < 1.0):
                    r += 0.01
                    continue
                found = True
                break
            except:
                pass
            r += 0.01
        
        if (found == False and last_good is None):
            return
        
        print('Fixed')
        new_expA, new_expB, new_expC = last_good
        new_splines = []
        found_spline = False
        self.expA = new_expA
        self.expB = new_expB
        self.expC = new_expC
        for spline in self.splines:
            if (r >=spline.r0 and r< spline.r1):
                b0 = spline.eval(r,0)
                b1 = spline.eval(r,1)
                b2 = spline.eval(r,2)/2.0
                b4 = spline.a4
                b3 = 4*b4*(r-spline.r0)+spline.a3
                new_splines.append(Spline4(r, spline.r1, b0, b1, b2, b3, b4))
                found_spline = True
            else:
                if (found_spline == False):
                    continue
                else:
                    new_splines.append(spline)
        self.splines = new_splines
        self.nsplines = len(self.splines)

    def __init__(self, name, knots, coeffs, calc_exp=True):
        self.name = name
        self.nsplines = len(knots)-1
        if len(coeffs) != self.nsplines:
            print(len(coeffs), self.nsplines)
            raise ValueError("Size inconsistent: knots and coeffs")
        if (len(coeffs[0]) == 4 or (self.nsplines == 1 and len(coeffs[0])==6)):
            self.order = 3
        elif (len(coeffs[0]) == 5):
            self.order = 4
        else:
            raise ValueError("Size inconsistent: number of coeffs")
        
        self.splines = []
        self.knots = knots
        self.cutoff = knots[-1]

        for i in range(0, len(knots)-1):
            r0 = knots[i]
            r1 = knots[i+1]
            ind = i
            if (self.order == 3):
                if (ind == len(coeffs)-1):
                    # l_coeffs = coeffs[ind]
                    # if (len(l_coeffs) <6):
                    #     l_coeffs += [0.0]*(6-len(l_coeffs))
                    # spline = Spline5(r0, r1, *l_coeffs)
                    spline = Spline5(r0, r1, *coeffs[ind])
                else:
                    spline = Spline3(r0, r1, *coeffs[ind])
            else:
                spline = Spline4(r0, r1, *coeffs[ind])
            self.splines.append(spline)
        
        
        if (calc_exp):
            try:
                self.expA = -2.0*self.splines[0].a2 / self.splines[0].a1
                self.expB = math.log(-self.splines[0].a1/self.expA) + self.expA*self.splines[0].r0
                self.expC = self.splines[0].a0 - math.exp(-self.expA*self.splines[0].r0+self.expB)
            except:
                print('Nag Expt.')
                self.expB = math.nan
                self.expC = math.nan
                self.fix_exponential()

    def __str__(self):
        with io.StringIO() as sout:
            if (self.order == 3):
                print('Spline', file=sout)
            else:
                print('Spline4', file=sout)
            print('{} {:.6f}'.format(self.nsplines, self.splines[-1].r1), file=sout)
            print('{:.6E} {:.6E} {:.6E}'.format(self.expA, self.expB, self.expC), file=sout)
            for spline in self.splines:
                print (str(spline), file=sout)

            return sout.getvalue()
    
    @classmethod
    def from_file(cls, filename):
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('Spline'):
                    line = next(f)
                    arr = line.split()
                    nspline = int(arr[0])
                    # cutoff = float(arr[1])
                    line = next(f)
                    arr = line.split()
                    exp_a = float(arr[0])
                    exp_b = float(arr[1])
                    exp_c = float(arr[2])
                    knots = []
                    coeffs = []
                    for _ in range(nspline):
                        line = next(f)
                        arr = line.split()
                        x0 = float(arr[0])
                        x1 = float(arr[1])
                        cofs = list(map(lambda x: float(x), arr[2:]))
                        coeffs.append(cofs)                    
                        if (len(knots) == 0):
                            knots.append(x0)    
                        knots.append(x1)
                    obj = cls(os.path.basename(filename), knots, coeffs, False)
                    obj.expA = exp_a
                    obj.expB = exp_b
                    obj.expC = exp_c
                    return obj
            raise RuntimeError('Not a valid repulsive potential.')


    

input_rep_1 = sys.argv[1]
input_rep_2 = sys.argv[2]

spl_rep_1 = RepulsivePotenial.from_file(input_rep_1)
spl_rep_2 = RepulsivePotenial.from_file(input_rep_2)

name = os.path.basename(input_rep_1).split('.')[0]
output_fig = '{}.eps'.format(name)
        
plt.title(name)
plt.ylim(-5, 5)

styles = ['c--', 'r--', 'g--']
labels = ['Rep. Potential', '1st derivative', '2nd derivative']
for i in range(3):
    xs = np.linspace(1.0, max(spl_rep_1.knots[-1], spl_rep_2.knots[-1]), num=400)
    ref_ys = []
    ys = []
    for r in xs:
        v = spl_rep_2.eval(r, i)
        v_ref = spl_rep_1.eval(r, i)
        ys.append(v)
        ref_ys.append(v_ref)
    
    if (i==0):
        error = abs(np.array(ref_ys)- np.array(ys))
        print( xs[np.argmax(error)], max(error) )


    plt.plot(xs, ref_ys, 'k-', linewidth=1)
    plt.plot(xs, ys, styles[i], linewidth=1, label=labels[i])
    
plt.legend(loc='upper right')
plt.xlabel('Diatomic Distance [a.u.]')
plt.ylabel('Energy [a.u.]')
plt.savefig(output_fig, format='eps')


    


                    