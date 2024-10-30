#!/usr/bin/env python3

from   scipy import *
from   scipy.optimize import leastsq

def inifit(x,y):
    p=[] ; a,b,c = polyfit(x,y,2)
    p.append(-b/(2*a))               # v0
    p.append(a*p[0]**2 + b*p[0] + c) # e0
    p.append((2*a*p[0]))             # b0
    p.append(2.)                     # bp
    return p

def bmeos(v,p):
    v0 = p[0] ; e0 = p[1] ; b0 = p[2] ; bp = p[3]  
    vv = (v0/v)**(2./3)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee    

def residuals(p,e,v):
    return e - bmeos(v,p)

def bmfit(x,y):
    p = [0,0,0,0] ; p = inifit(x,y)
    return leastsq(residuals, p, args=(y, x))


def sortstrain(s,e):
    ss=[] ; ee=[] ; ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee

def eos(v,v0,e0,b0,bp):
    vv = (v0/v)**(2./3)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee   

def calculate_eos(volume, energy):
    factor = 2.0
    # bohr_radius     = 0.529177
    # joule2hartree   = 4.3597482
    # unitconv        = joule2hartree/bohr_radius**3*10.**3*factor

    electron_charge = 1.602176565e-19 
    bohr_radius     = 5.2917721092e-11
    rydberg2ev      = 13.605698066      
    unitconv        = electron_charge*rydberg2ev/(1e9*bohr_radius**3)*factor

    volume, energy = sortstrain(volume,energy)

    p = bmfit(volume,energy)

    v0 = p[0][0] ; e0 = p[0][1] ; b0 = p[0][2]*unitconv ; bp =  p[0][3]

    

    chi = 0 ; ebm = [] ; fchi = 0
    for i in range(len(volume)): 
        chi=chi+residuals(p[0],energy[i],volume[i])**2
        fchi=fchi+(energy[i]-eos(volume[i],v0,e0,b0/unitconv,bp))**2
        ebm.append(eos(volume[i],v0,e0,b0/unitconv,bp))
        

    chi=sqrt(chi)/len(volume)

    res = {}
    res['V0'] = v0
    res['B0'] = b0
    res['BP'] = bp
    res['chi'] = log10(chi)

    return res