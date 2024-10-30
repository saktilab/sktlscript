#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   sys   import stdin
from   math  import sqrt
from   math  import factorial
import matplotlib
matplotlib.use('PS')
from   matplotlib.pylab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk 
import matplotlib.pylab  as pyl
import numpy
import sys
import os  
import math 
import numpy
from   scipy import *
from   scipy.optimize import leastsq

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

def printderiv(damax,bb,nfit,output_file):
    print ('%16.10f'%(damax), file=output_file)
    for j in range(0,nfit):
        if (bb[j] != 99999999):
            print ('%14.6f'%(bb[j]), end=" ", file=output_file)
    print ( file=output_file)
    return 

def printnice(etam,bb,nfit,nderiv,orderstep,strainpoints,dc,dl):
    print ()
    print ("###########################################\n")
    print ("Fit data-----------------------------------\n")
    print ("Deformation code             ==>", dc)  ##'%2i'%(dc))
    print ("Deformation label            ==>", dl)
    print ("Maximum value of the strain  ==>", '%10.8f'%(etam)) 
    print ("Number of strain values used ==>", '%2i'%(strainpoints),"\n") 
    print ("Fit results for the derivative of order", '%3i'%(nderiv),"\n")  
    for j in range(0,nfit):
        order=nderiv+j*orderstep
        if (bb[j] != 99999999):
            print ("Polynomial of order", '%2i'%(order), "==>", end=" ") 
            print ('%10.2f'%(bb[j]), "[GPa]")
    print ()
    return 

def sortstrain(s,e):
    ss=[] ; ee=[] ; ww=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee
    
#-------------------------------------------------------------------------------

def shell_value(variable,vlist,default):
    v = default
    e = False
    for i in range(len(vlist)):
        if ( vlist[i] == variable ):
            v = os.environ[variable] 
            e = True 
            break
    return v, e
    
def eos(v,v0,e0,b0,bp):
    vv = (v0/v)**(2./3)
    ff = vv - 1.
    ee = e0 + 9.*b0*v0/16. * ( ff**3*bp + ( 6.-4.*vv )*ff**2 )
    return ee   

#-------------------------------------------------------------------------------

factor     = 2.
startorder = 0

lrydberg = os.path.exists('quantum-espresso') or os.path.exists('vasp')
lplanar  = os.path.exists('planar')
lstartor = os.path.exists('startorder')
lexsinfo = os.path.exists('INFO-elastic-constants')
lexsev   = os.path.exists('energy-vs-volume')

if (lrydberg): factor=1.0

if (lplanar):   
    input_planar = open('planar',"r")
    factor=factor*float(input_planar.readline().strip().split()[0])

if (lstartor): 
    input_startorder = open('startorder',"r")
    startorder=int(input_startorder.readline().strip().split()[0])

if (not(lexsev)):   sys.exit("ERROR: file energy-vs-alat not found!\n")

#-------------------------------------------------------------------------------

bohr_radius     = 0.529177
joule2hartree   = 4.3597482
joule2rydberg   = joule2hartree/2.
unitconv        = joule2hartree/bohr_radius**3*10.**3*factor

electron_charge = 1.602176565e-19 
bohr_radius     = 5.2917721092e-11
rydberg2ev      = 13.605698066      
unitconv        = electron_charge*rydberg2ev/(1e9*bohr_radius**3)*factor

#-------------------------------------------------------------------------------

input_energy = open('energy-vs-volume',"r")
output_birch = open('birch-murnaghan',"w")

energy = [] ; volume = []

while True:
    line = input_energy.readline().strip()
    if len(line) == 0: break
    energy.append(float(line.split()[1]))
    volume.append(float(line.split()[0]))

nv = len(volume)
if (nv < 4): sys.exit("\nERROR: Too few volumes ("+str(nv)+")!\n")

#-------------------------------------------------------------------------------

volume, energy = sortstrain(volume,energy)

#-------------------------------------------------------------------------------

p = bmfit(volume,energy)

v0 = p[0][0] ; e0 = p[0][1] ; b0 = p[0][2]*unitconv ; bp =  p[0][3]

a0sc=(1*p[0][0])**(0.33333333333)
abcc=(2*p[0][0])**(0.33333333333)
afcc=(4*p[0][0])**(0.33333333333)
ahcp=(p[0][0]/math.cos(math.pi/6)/1.633)**(0.3333333333333)

chi = 0 ; ebm = [] ; eee = [] ; fchi = 0
for i in range(len(volume)): 
    chi=chi+residuals(p[0],energy[i],volume[i])**2
    fchi=fchi+(energy[i]-eos(volume[i],v0,e0,b0/unitconv,bp))**2
    ebm.append(eos(volume[i],v0,e0,b0/unitconv,bp))
    #eee.append(eos(volume[i],v0,e0,378./unitconv,bp))
    #print volume[i], ebm[i]

chi=sqrt(chi)/len(volume)

fmt='%10.4f'
amt='%10.4f'
bmt='%8.3f'
pmt='%16.10f'
lmt='%10.2f'

print ()
print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print ("                                          Atomic units")
print ("     V0        B0        BP         a-sc       a-bcc      a-fcc      a-hcp     log(chi)")
print (fmt%(v0), bmt%(b0), bmt%(bp)," ",end="")  
print (amt%(a0sc), amt%(abcc), amt%(afcc), amt%(ahcp), lmt%(log10(chi)))
print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print ()
print (pmt%(v0), pmt%(b0), pmt%(bp),end="", file=output_birch)  
print (pmt%(a0sc), pmt%(abcc), pmt%(afcc), pmt%(ahcp),file=output_birch)
print (pmt%(log10(chi)),file=output_birch)

print ()
print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print ("                                          Angstrom")
print ("     V0        B0        BP         a-sc       a-bcc      a-fcc      a-hcp     log(chi)")
print (fmt%(v0*0.529177**3), bmt%(b0), bmt%(bp)," ",end="")  
print (amt%(a0sc*0.529177), amt%(abcc*0.529177), amt%(afcc*0.529177), amt%(ahcp*0.529177), lmt%(log10(chi)))
print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print ()
print (fmt%(v0*0.529177**3), bmt%(b0), bmt%(bp)," ",end="", file=output_birch)  
print (pmt%(a0sc*0.529177), pmt%(abcc*0.529177), pmt%(afcc*0.529177), pmt%(ahcp*0.529177), lmt%(log10(chi)),file=output_birch)

input_energy.close()
output_birch.close()

#*******************************************************************************

current = os.environ['PWD']
ev_list = os.environ.keys()

# rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
# rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
# showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
# dpipng = int(shell_value('DPIPNG',ev_list,300)[0])
dpipng = 300

#-------------------------------------------------------------------------------

xlabel = u'Volume [Bohr\u00B3]'
ylabel = r'Energy [Ha]'
if (os.path.exists('quantum-espresso')): ylabel = r'Energy [Ry]'
if (os.path.exists('vasp')): ylabel = r'Energy [Ry]'

#-------------------------------------------------------------------------------

xvol = numpy.linspace(volume[0],volume[-1],100)

xene = [] 
for i in range(len(xvol)):
    xene.append(eos(xvol[i],v0,e0,b0/unitconv,bp))    

#-------------------------------------------------------------------------------

fontlabel=20
fonttick=16

params = {'ytick.minor.size': 6,
          'xtick.major.pad': 8,
          'ytick.major.pad': 4,
          'patch.linewidth': 2.,
          'axes.linewidth': 2.,
          'lines.linewidth': 1.8,
          'lines.markersize': 8.0,
          'axes.formatter.limits': (-4, 6)}

plt.rcParams.update(params)
plt.subplots_adjust(left=0.21, right=0.93,
                    bottom=0.18, top=0.88,
                    wspace=None, hspace=None)
                           
yfmt = ptk.ScalarFormatter(useOffset=True,useMathText=True)                           
                           
figure = plt.figure(1, figsize=(8,5.5))  
ax     = figure.add_subplot(111)
ax.text(0.5,-0.13,xlabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center')
ax.text(-0.19,0.5,ylabel,size=fontlabel,
        transform=ax.transAxes,ha='center',va='center',rotation=90)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)
plt.xticks(size=fonttick)
plt.yticks(size=fonttick)
pyl.grid(True)

plt.plot(xvol,xene,'b-',label='birch-murnaghan fit')
plt.plot(volume,energy,'go',label='calculated')
plt.plot(v0,e0,'ro')
print(v0, e0)
plt.legend(loc=9,borderaxespad=.8,numpoints=1)

ymax  = max(max(xene),max(energy))
ymin  = min(min(xene),min(energy))
dxx   = abs(max(volume)-min(volume))/18
dyy   = abs(ymax-ymin)/18
ax.yaxis.set_major_formatter(yfmt)
ax.set_xlim(min(volume)-dxx,max(volume)+dxx)
ax.set_ylim(ymin-dyy,ymax+dyy)

ax.xaxis.set_major_locator(MaxNLocator(7))

plt.savefig('PLOT.ps', orientation='portrait',format='eps')
plt.savefig('PLOT.png',orientation='portrait',format='png',dpi=dpipng)

#-------------------------------------------------------------------------------





