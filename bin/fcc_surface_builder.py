#!/usr/bin/env python
# Purpose: Calculate FCC (111), (110), and (100) surface energies
# Author: Sungho Kim and Laalitha Liyanage

import os
import sys
import math

usage="""
        Usage: ./gen_fcc_surface.py a surf vacuum nx ny nz adatom

        Mandatory arguments
        -------------------
        a - equilibrium lattice constant
        surf - Type of surface 100, 110 or 111

        Optional arguments
        -------------------
        vacuum - length of vacuum; DEFAULT = 8.0 angstroms
        nx,ny,nz - periodicity of supercell; DEFAULT (1,1,1)
        adatom - 1/0 (True/False); DEFAULT = 0 (False)
"""


#Default setting
#--------------------------------------------------------------------
vacuum = 8.0
adatom = 0

#--------------------------Surface (100)-----------------------------

def gen_data_for_fcc(a,nx=2,ny=2,nz=4):
  """ Generate datafile of FCC structure with lattice constant a """
  xa=[]; ya=[]; za=[]
  x0 = 0.0
  bx,by,bz = a*nx,a*ny,a*nz+vacuum
  x,y,z = bx,by,bz

  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        xa.append(0 + i*a); ya.append(  0 + j*a); za.append(  0 + k*a)
        xa.append(  0 + i*a); ya.append(a/2 + j*a); za.append(a/2 + k*a)
        xa.append(a/2 + i*a); ya.append(  0 + j*a); za.append(a/2 + k*a)
        xa.append(a/2 + i*a); ya.append(a/2 + j*a); za.append(  0 + k*a)
  if adatom != 0:
#        xa.append(x0); ya.append(x0); za.append(x0+nz*a)
        xa.append(bx/2.); ya.append(by/2.); za.append(x0+nz*a)
  return xa,ya,za,bx,by,bz

#--------------------------Surface (110)-----------------------------

def gen_data_for_110_fcc(a,nx=4,ny=2,nz=1):
  """ Generate datafile of FCC surface: 110:x, 112:y, 111:z """
  xa=[]; ya=[]; za=[]
  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4. * a
  y2 = math.sqrt(6)/4. * a
  y3 = math.sqrt(6)/6. * a
  y4 = math.sqrt(6)*5./12. * a
  y5 = math.sqrt(6)*2./6. * a
  y6 = math.sqrt(6)/12 * a
  z3 = math.sqrt(3)/3. * a
  z5 = math.sqrt(3)*2./3. * a
  bx,by,bz = ax*nx + vacuum, ay*ny, az*nz
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(x0+k*az)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(x0+k*az)
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(z3+k*az)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(z3+k*az)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(z5+k*az)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(z5+k*az)
  if adatom != 0:
        xa.append(x0+nx*ax); ya.append(by/2.); za.append(bz/2.)
  return xa,ya,za,bx,by,bz

#--------------------------Surface (111)-----------------------------

def gen_data_for_111_fcc(a,nx=2,ny=2,nz=4):
  """ Generate datafile of FCC surface: 110:x, 112:y, 111:z """
  xa=[]; ya=[]; za=[]
  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4 * a
  y2 = math.sqrt(6)/4 * a
  y3 = math.sqrt(6)/6 * a
  y4 = math.sqrt(6)*5/12 * a
  y5 = math.sqrt(6)*2/6 * a
  y6 = math.sqrt(6)/12 * a
  bx,by,bz = ax*nx, ay*ny, az*nz+vacuum
  for i in range(nx):
    for j in range(ny):
      layer = 0
      for k in range(nz):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1
  if adatom != 0:
        xa.append(bx/2.); ya.append(by/2.); za.append(x0+nz*az)
  return xa,ya,za,bx,by,bz
#----------------------------POSCAR generation------------------------------------------------
def gen_poscar(xa,ya,za,bx,by,bz):
  fout = open("POSCAR","w")
  fout.write("Fe\n")
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(bx,0,0))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(0,by,0))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(0,0,bz))
  fout.write("%d\n"%len(xa))
#  fout.write("Selective Dynamics\n")
  fout.write("Cart\n")
  for i in range(len(xa)):
    fout.write("%22.16f %22.16f %22.16f\n"%(xa[i],ya[i],za[i]))
#    fout.write("%22.16f %22.16f %22.16f F F T\n"%(xa[i],ya[i],za[i]))
  fout.close()
  return len(xa)

#-------------------------------Main program---------------------------------------------------

if len(sys.argv) > 2:
        if len(sys.argv) == 3:
                a_latt = float(sys.argv[1])
                surf = sys.argv[2]
                if surf == '100' :
                        xa,ya,za,bx,by,bz = gen_data_for_fcc(a_latt)
                        gen_poscar(xa,ya,za,bx,by,bz)

                elif surf == '110':
                        xa,ya,za,bx,by,bz = gen_data_for_110_fcc(a_latt)
                        gen_poscar(xa,ya,za,bx,by,bz)

                elif surf == '111':
                        xa,ya,za,bx,by,bz = gen_data_for_111_fcc(a_latt)
                        gen_poscar(xa,ya,za,bx,by,bz)

        elif len(sys.argv) == 8:
                a_latt = float(sys.argv[1])
                surf = sys.argv[2]
                vacuum = float(sys.argv[3])
                nx = int(sys.argv[4])
                ny = int(sys.argv[5])
                nz = int(sys.argv[6])
                adatom = int(sys.argv[7])

                if surf == '100' :
                        xa,ya,za,bx,by,bz = gen_data_for_fcc(a_latt,nx,ny,nz)
                        gen_poscar(xa,ya,za,bx,by,bz)

                elif surf == '110':
                        xa,ya,za,bx,by,bz = gen_data_for_110_fcc(a_latt,nx,ny,nz)
                        gen_poscar(xa,ya,za,bx,by,bz)

                elif surf == '111':
                        xa,ya,za,bx,by,bz = gen_data_for_111_fcc(a_latt,nx,ny,nz)
                        gen_poscar(xa,ya,za,bx,by,bz)



else:
    print ("Error: wrong number of arguments!!!")
    print (usage)
