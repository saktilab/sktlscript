#!/usr/bin/env /Users/adit/miniconda3/bin/python

import sys
import argparse
from ase.units import Ha, Bohr
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.offsite_twocenter import Offsite2cTable
from ase.data import atomic_numbers, covalent_radii

parser = argparse.ArgumentParser(description='Program to calculate the atomic eigen value and Hubbard value')
parser.add_argument('-e','--element',type=str,required=True,help='Type of element you want to calculate')
parser.add_argument('-f','--functional',type=str,default='LDA',help='Type of exchange functional. Options: LDA, GGA')
parser.add_argument('-c','--conf',type=str,required=True,help='Electronic configuration of the element. Example: [Ne] 3s2 3p2')
parser.add_argument('-v','--valence',type=str,required=True, help="Valence shell of the element. Example: 3s,3p")
parser.add_argument('-r','--rel',type=str,default='False',help='Whether the scalar relativistic calculations will be performed. Default: False.')
parser.add_argument('-p','--perturb',type=str,default='False',help='Whether the perturbative confinement will be included. Default: False.')
parser.add_argument('-r0',type=str,default='None',help='r0 value of the power confinement. Default: 40.0')
parser.add_argument('-s',type=str,default='None',help='s value of the power confinement. Default: 4')
parser.add_argument('-super',type=str,default='density',help='Type of superposition in electronic parameterization. Default: density.')
parser.add_argument('-maxiter',type=int,default=500,help='Maximum number of charge iterations. Default: 500.')
parser.add_argument('-mix',type=float,default=0.05,help='Mixing coefficients for the charge iterations. Default: 0.05.')

opt = parser.parse_args(sys.argv[1:])
elements = opt.element.split(",")
# Extracting electronic configurations and valences for each element
conf = opt.conf.split(",")
val_ele = opt.valence.split(";")
if val_ele[0] == '1s':
    val_ele0 = [val_ele[0]]
else:
    val_ele0 = val_ele[0].split(",")
if len(val_ele) != 1:
    if val_ele[1] == '1s':
        val_ele1 = [val_ele[1]]
    else:
        val_ele1 = val_ele[1].split(",")
# valence = dict(zip(elements,val_ele))
# elec_conf = dict(zip(elements,conf))

if len(elements) !=1:
#     atom1 = AtomicDFT('{}'.format(elements[0]), 
#                  xc='{}'.format(opt.functional), 
#                  configuration='{}'.format(conf[0]),
#                  valence=val_ele0,
#                  scalarrel='{}'.format(opt.rel),
#                  perturbative_confinement='{}'.format(opt.perturb),
#                  txt='-'
# )
# # Run the calculations
#     atom1.run()
#     atom1.plot_Rnl('{}_Rnl_free.png'.format(elements[0])) # Plot the radial parts of the valence orbitals
#     atom1.plot_rho('{}_rho_free.png'.format(elements[0])) # Plot the valence orbital densities and total electron density
#     print(10*'=')
#     atom2 = AtomicDFT('{}'.format(elements[1]), 
#                  xc='{}'.format(opt.functional), 
#                  configuration='{}'.format(conf[1]),
#                  valence=val_ele1,
#                  scalarrel='{}'.format(opt.rel),
#                  perturbative_confinement='{}'.format(opt.perturb),
#                  txt='-'
# )
# # Run the calculations
#     atom2.run()
#     atom2.plot_Rnl('{}_Rnl_free.png'.format(elements[1])) # Plot the radial parts of the valence orbitals
#     atom2.plot_rho('{}_rho_free.png'.format(elements[1])) # Plot the valence orbital densities and total electron density
#     print(10*'=')

#     print("Details for {}".format(elements[0]))
#     for nl in atom1.valence:
#         e = atom1.get_eigenvalue(nl)
#         print(nl, e, 'Hartree =', e*Ha, 'eV')

#     print("Details for {}".format(elements[1]))
#     for nl in atom2.valence:
#         e = atom2.get_eigenvalue(nl)
#         print(nl, e, 'Hartree =', e*Ha, 'eV')
    if opt.r0 != 'None':
        atom1 = AtomicDFT('{}'.format(elements[0]), 
                 xc='{}'.format(opt.functional), 
                 configuration='{}'.format(conf[0]),
                 valence=val_ele0,
                 scalarrel='{}'.format(opt.rel),
                 maxiter=opt.maxiter,
                 mix=opt.mix,
                 perturbative_confinement='{}'.format(opt.perturb),
                 confinement=PowerConfinement(r0=float(opt.r0), s=float(opt.s)),
                 txt=None
    )
        atom1.run()
        atom1.plot_Rnl('{}_Rnl_conf.png'.format(elements[0]))
        atom1.plot_rho('{}_rho_conf.png'.format(elements[0]))
        
        atom2 = AtomicDFT('{}'.format(elements[1]), 
                 xc='{}'.format(opt.functional), 
                 configuration='{}'.format(conf[1]),
                 valence=val_ele1,
                 maxiter=opt.maxiter,
                 mix=opt.mix,
                 scalarrel='{}'.format(opt.rel),
                 perturbative_confinement='{}'.format(opt.perturb),
                 confinement=PowerConfinement(r0=float(opt.r0), s=float(opt.s)),
                 txt=None
    )
        atom2.run()
        atom2.plot_Rnl('{}_Rnl_conf.png'.format(elements[1]))
        atom2.plot_rho('{}_rho_conf.png'.format(elements[1]))        
        rmin, dr, N = 0.4, 0.02, 600
        off2c = Offsite2cTable(atom1,atom2)
        off2c.run(rmin,dr,N,superposition=opt.super,xc=opt.functional)
        off2c.write()
        off2c.plot('{}-{}_offsite2c.png'.format(elements[0],elements[1]))
        print("Done!")        
    
        
    else:
        rcov0 = covalent_radii[atomic_numbers[elements[0]]]/Bohr
        rcov1 = covalent_radii[atomic_numbers[elements[1]]]/Bohr
        confine0 = PowerConfinement(r0=3*rcov0,s=2)
        confine1 = PowerConfinement(r0=3*rcov1,s=2)
        pow1 = []
        pow2 = []
        for i in val_ele0:
            pow1.append(PowerConfinement(r0=2*rcov0,s=2))
        for i in val_ele1:
            pow2.append(PowerConfinement(r0=2*rcov1,s=2))
        wf_conf1 = dict(zip(val_ele0,pow1))
        wf_conf2 = dict(zip(val_ele1,pow2))
        atom1 = AtomicDFT('{}'.format(elements[0]), 
                 xc='{}'.format(opt.functional), 
                 configuration='{}'.format(conf[0]),
                 valence=val_ele0,
                 maxiter=opt.maxiter,
                 mix=opt.mix,
                 scalarrel='{}'.format(opt.rel),
                 confinement=confine0,
                 wf_confinement=wf_conf1,
                 perturbative_confinement='{}'.format(opt.perturb),
                 txt=None
    )
        atom1.run()
        atom1.plot_Rnl('{}_Rnl_conf.png'.format(opt.element))
        atom1.plot_rho('{}_rho_conf.png'.format(opt.element))

        atom2 = AtomicDFT('{}'.format(elements[1]), 
                 xc='{}'.format(opt.functional), 
                 configuration='{}'.format(conf[1]),
                 valence=val_ele1,
                 scalarrel='{}'.format(opt.rel),
                 maxiter=opt.maxiter,
                 mix=opt.mix,
                 confinement=confine1,
                 wf_confinement=wf_conf2,
                 perturbative_confinement='{}'.format(opt.perturb),
                 txt=None
    )
        atom2.run()
        atom2.plot_Rnl('{}_Rnl_conf.png'.format(opt.element))
        atom2.plot_rho('{}_rho_conf.png'.format(opt.element))        
        rmin, dr, N = 0.4, 0.02, 600
        off2c = Offsite2cTable(atom1,atom2)
        off2c.run(rmin,dr,N,superposition=opt.super,xc=opt.functional)
        off2c.write()
        off2c.plot('{}-{}_offsite2c.png'.format(elements[0],elements[1]))
        print("Done!")    
else:
#     atom1 = AtomicDFT('{}'.format(elements[0]), 
#                  xc='{}'.format(opt.functional), 
#                  configuration='{}'.format(conf[0]),
#                  valence=val_ele0,
#                  scalarrel='{}'.format(opt.rel),
#                  perturbative_confinement='{}'.format(opt.perturb),
#                  txt='-'
# )
# # Run the calculations
#     atom1.run()
#     atom1.plot_Rnl('{}_Rnl_free.png'.format(elements[0])) # Plot the radial parts of the valence orbitals
#     atom1.plot_rho('{}_rho_free.png'.format(elements[0])) # Plot the valence orbital densities and total electron density
#     print(10*'=')
    
    # print("Details for {}".format(elements[0]))
    # for nl in atom1.valence:
    #     e = atom1.get_eigenvalue(nl)
    #     print(nl, e, 'Hartree =', e*Ha, 'eV')

    if opt.r0 != 'None':
        atom1 = AtomicDFT('{}'.format(elements[0]), 
                 xc='{}'.format(opt.functional), 
                 configuration='{}'.format(conf[0]),
                 valence=val_ele0,
                 scalarrel='{}'.format(opt.rel),
                 maxiter=opt.maxiter,
                 mix=opt.mix,
                 perturbative_confinement='{}'.format(opt.perturb),
                 confinement=PowerConfinement(r0=float(opt.r0),s=float(opt.s)),
                 txt=None
    )
        atom1.run()
        atom1.plot_Rnl('{}_Rnl_conf.png'.format(elements[0]))
        atom1.plot_rho('{}_rho_conf.png'.format(elements[0]))
        rmin, dr, N = 0.4, 0.02, 600
        off2c = Offsite2cTable(atom1,atom1)
        off2c.run(rmin,dr,N,superposition=opt.super,xc=opt.functional)
        off2c.write()
        off2c.plot('{}-{}_offsite2c.png'.format(elements[0],elements[0]))
        print("Done!")    
        
    else:
        rcov = covalent_radii[atomic_numbers[opt.element]]/Bohr
        confine = PowerConfinement(r0=3*rcov,s=2)
        pow1 = []
        for i in val_ele0:
            pow1.append(PowerConfinement(r0=2*rcov,s=2))

        wf_conf1 = dict(zip(val_ele0,pow1))
        # wf_conf2 = dict(zip(val_ele1,values))
        atom1 = AtomicDFT('{}'.format(elements[0]), 
                 xc='{}'.format(opt.functional), 
                 configuration='{}'.format(conf[0]),
                 valence=val_ele0,
                 maxiter=opt.maxiter,
                 mix=opt.mix,
                 scalarrel='{}'.format(opt.rel),
                 confinement=confine,
                 wf_confinement=wf_conf1,
                 perturbative_confinement='{}'.format(opt.perturb),
                 txt=None
    )
        atom1.run()
        atom1.plot_Rnl('{}_Rnl_conf.png'.format(opt.element))
        atom1.plot_rho('{}_rho_conf.png'.format(opt.element)) 
        rmin, dr, N = 0.4, 0.02, 600
        off2c = Offsite2cTable(atom1,atom1)
        off2c.run(rmin,dr,N,superposition=opt.super,xc=opt.functional)
        off2c.write()
        off2c.plot('{}-{}_offsite2c.png'.format(elements[0],elements[0]))
        print("Done!")
            