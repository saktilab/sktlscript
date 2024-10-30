#!/usr/bin/env python3

import sys
import math
import io
import decimal
import os
import numpy as np
import scipy.interpolate

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

    def __init__(self, name, knots, coeffs):
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

        for i in range(0, len(knots)-1):
            r0 = knots[i]
            r1 = knots[i+1]
            ind = i
            if (self.order == 3):
                if (ind == len(coeffs)-1):
                    spline = Spline5(r0, r1, *coeffs[ind])
                else:
                    spline = Spline3(r0, r1, *coeffs[ind])
            else:
                spline = Spline4(r0, r1, *coeffs[ind])
            self.splines.append(spline)
        
        

        self.expA = -2.0*self.splines[0].a2 / self.splines[0].a1
        self.expB = math.log(-self.splines[0].a1/self.expA) + self.expA*self.splines[0].r0
        self.expC = self.splines[0].a0 - math.exp(-self.expA*self.splines[0].r0+self.expB)

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
                if line.startswith('Spline4'):
                    line = next(f)
                    arr = line.split()
                    nspline = int(arr[0])
                    cutoff = float(arr[1])
                    line = next(f)
                    arr = line.split()
                    exp_a = float(arr[0])
                    exp_b = float(arr[1])
                    exp_c = float(arr[2])
                    knots = []
                    coeffs = []
                    for i in range(nspline):
                        line = next(f)
                        arr = line.split()
                        x0 = float(arr[0])
                        x1 = float(arr[1])
                        cofs = list(map(lambda x: float(x), arr[2:]))
                        coeffs.append(cofs)                    
                        if (len(knots) == 0):
                            knots.append(x0)    
                        knots.append(x1)
                    return cls(os.path.basename(filename), knots, coeffs)
            raise RuntimeError('Not a valid repulsive potential.')

if __name__ == '__main__':
    ori_rep_file = sys.argv[1]
    starting_dis = float(sys.argv[2])
    shift = float(sys.argv[3]) ##shift in hartree
    interval = 0.02

    rep = RepulsivePotenial.from_file(ori_rep_file)
    start = rep.knots[0]
    end = rep.knots[-1]
    

    
    numpt = int((end-start)/interval)+1
    new_knots = np.linspace(start, end, numpt)
    new_ys = []
    new_xs = [] 
    final_coeffs = None
    final_x0 = None 
    for x_val in new_knots:
        new_y = rep.eval(x_val)
        new_y += shift
        if (x_val >= starting_dis):
            
            new_xs.append(x_val)
            new_ys.append(new_y)

            final_x0 = x_val
            start_2nd = rep.eval(new_xs[0], 2)
            end_2nd = rep.eval(new_xs[-1], 1)*0.93
            
            res_spline = scipy.interpolate.CubicSpline(new_xs, new_ys, bc_type=((2, start_2nd), (1, end_2nd)))

            # print(res_spline(final_x0, 0), res_spline(final_x0, 1), res_spline(final_x0, 2))
           
            dx = end - x_val
            alpha = res_spline(x_val) 
            beta = res_spline(x_val, 1)
            gamma = res_spline(x_val, 2)

            mat = np.zeros((6,6))
            mat[0][0] = 1.0
            mat[0][1] = dx
            mat[0][2] = dx**2.0
            mat[0][3] = dx**3.0
            mat[0][4] = dx**4.0
            mat[0][5] = dx**5.0

            mat[1][0] = 0.0
            mat[1][1] = 1.0
            mat[1][2] = 2.0*dx
            mat[1][3] = 3.0*dx**2.0
            mat[1][4] = 4.0*dx**3.0
            mat[1][5] = 5.0*dx**4.0

            mat[2][0] = 0.0
            mat[2][1] = 0.0
            mat[2][2] = 2.0
            mat[2][3] = 6.0*dx
            mat[2][4] = 12.0*dx**2.0
            mat[2][5] = 20.0*dx**3.0

            mat[3][0] = 1.0
            mat[4][1] = 1.0
            mat[5][2] = 2.0

            v = np.zeros(6)
            v[3] = alpha
            v[4] = beta
            v[5] = gamma

            final_coeffs = np.linalg.solve(mat, v)
            # print(final_coeffs)

            solved_coeffs = [[x[3], x[2], x[1], x[0]] for x in np.transpose(res_spline.c).tolist()]
            solved_coeffs.append(final_coeffs)
            new_xs.append(end)
            
            # print(new_xs)
            # print(solved_coeffs, len(solved_coeffs))
            # print(len(new_xs))
            final_rep = RepulsivePotenial(rep.name, new_xs, solved_coeffs)
            print(final_rep)
            # print(final_coeffs[0] + dx*final_coeffs[1]+final_coeffs[2]*dx*dx+final_coeffs[3]*dx*dx*dx+final_coeffs[4]*dx*dx*dx*dx+final_coeffs[5]*dx*dx*dx*dx*dx)
            # print('{:.5f} {:20.12f} {:20.12f} {:20.12f}'.format(x_val, new_y, new_y2, damping))
            break
        else:
            new_xs.append(x_val)
            new_ys.append(new_y)
    
    


            


        
        

    


                    