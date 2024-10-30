#!/usr/bin/env python3

import sys
import argparse
import pathlib
import subprocess
import os
import shutil
import gmxscript
import numpy as np

parser = argparse.ArgumentParser(description='Waterbox generator')
parser.add_argument('-n', '--numWater', default=128, type=int, help='Number of water')
parser.add_argument('-o', '--output', type=str, default='waterbox.dcdftb')
parser.add_argument('-s', '--subsystem', action='store_true', help='generate subsystem')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-t', '--temperature', default=300.0, type=float, help='Temperature(K)')
group.add_argument('-l', '--boxsize', type=float, help='Box size(Angstrom)')
group.add_argument('-d', '--density', type=float, help='Density(g/cm^3)')

parser.add_argument('--packmol', type=str, default='/Users/radenwibawasakti/opt/packmol/exe/packmol')

opt = parser.parse_args(sys.argv[1:])

#T in K
#Formula from https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
def water_density(T):
    vK = T-273.15 
    return -9.204453627e-11*vK**4+3.420742008672e-8*vK**3-7.08919807166417e-6*vK*vK+4.375294545181970e-5*vK+0.99988826440573500000

def run_command(command, args=[], input_str=None):
    command_list = [command]
    if (len(args) > 0):
        command_list.extend(args)

    if (input_str is not None):
        res =  subprocess.run(command_list, stdout=subprocess.PIPE, input=input_str.encode(encoding='ascii'))
    else:
        res =  subprocess.run(command_list, stdout=subprocess.PIPE)
    return res.stdout.decode()


water_mass = 18.01528  #(amu/mol)
avogadro = 6.02214086e+23 # mol^-1


single_water = """ATOM      1  OW  SOL   243      16.510  17.500  12.870  1.00  0.00
ATOM      2  HW1 SOL   243      16.590  17.680  11.940  1.00  0.00
ATOM      3  HW2 SOL   243      15.910  18.160  13.200  1.00  0.00
ATOM      4  MW  SOL   243      16.460  17.570  12.820  1.00  0.00
"""

tip4p_fb_itp = """; TIP4P-FB water model:
; Lee-Ping Wang, Todd J. Martinez and Vijay S. Pande. Building force fields - an automatic, systematic and reproducible approach.  
; Journal of Physical Chemistry Letters, 2014, 5, pp 1885-1891.  DOI:10.1021/jz500737m


[ defaults ]
1    2    yes    0.5    0.8333

[ atomtypes ]
OW_tip4pfb    8     15.99940     0.00000     A    3.16555e-01  7.49279e-01
HW_tip4pfb    1      1.00800     0.00000     A    0.00000e+00  0.00000e+00
MW_tip4pfb    0      0.00000     0.00000     D    0.00000e+00  0.00000e+00 ; Same as other virtual sites

[ moleculetype ]
SOL        2

[ atoms ]
     1  OW_tip4pfb   1    SOL     OW      1       0.00000
     2  HW_tip4pfb   1    SOL    HW1      1       0.52587
     3  HW_tip4pfb   1    SOL    HW2      1       0.52587
     4  MW_tip4pfb   1    SOL     MW      1      -1.05174

#ifndef FLEXIBLE

[ settles ]
1    1    0.09572    0.15139

#else

[ bonds ]
; Copied straight from amber99sb-ildn.ff/tip4pew.itp.
; This is a rigid water model - do NOT use flexible parameters
1    2    1    0.09572    502416.0
1    3    1    0.09572    502416.0
        
[ angles ]
; Copied straight from amber99sb-ildn.ff/tip4pew.itp.
; This is a rigid water model - do NOT use flexible parameters
2    1    3    1    104.52   628.02

#endif

[ virtual_sites3 ]
4    1    2    3    1       0.0898426712735     0.0898426712735

[ exclusions ]
1    2    3    4
2    1    3    4
3    1    2    4
4    1    2    3

[ system ]
Built with Packmol

[ molecules ]
SOL               {}
"""

minim_mdp = """
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 100.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = {}       ; Short-range electrostatic cut-off
rvdw            = {}       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
;define = -DFLEXIBLE
"""

nvt_mdp="""
define                  = -DFLEXIBLE
integrator              = md        ; leap-frog integrator
nsteps                  = 50000    ; 1 * 500000 = 50 ps
dt                      = 0.001     ; 1 fs
nstxout                 = 100       ; save coordinates every 1.0 ps
;nstenergy               = 100       ; save energies every 1.0 ps
;nstlog                  = 100       ; update log file every 1.0 ps

; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = {}       ; short-range electrostatic cutoff (in nm)
rvdw                    = {}       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = System   ; two coupling groups - more accurate
tau_t                   = 0.1               ; time constant, in ps
ref_t                   = 300              ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
"""

n_water = opt.numWater
mol_water = n_water/avogadro
temperature = opt.temperature

if opt.temperature:
    density = water_density(temperature)
if opt.density:
    density = opt.density

box_size = (water_mass*n_water*1.0E24/(density*avogadro))**(1.0/3.0)

if (opt.boxsize):
    box_size = opt.boxsize
    density = (water_mass*n_water*1.0E24)/(box_size**3.0*avogadro)

if opt.temperature:
    print (f'Temperature: {temperature:12.4f} K')
print     (f'Density:     {density:12.6f} g/cm^3')
print     (f'Box length:  {box_size:12.4f} AA')



if box_size < 12:
    print('The box is too small')
    sys.exit(1)

half_box = min(1.0, box_size/20.0-0.1)
print(half_box)

def run_prep_box():
    cwd = os.path.abspath(os.curdir)
    print('Making pure solvent system using packmol')
    os.makedirs('prep', exist_ok=True)
    os.chdir('prep')

    with open('water.pdb', 'w') as f:
        print(single_water, file=f)

    with open('packmol.inp', 'w') as fout:
        print("""
        tolerance 2.5
        filetype pdb
        output system_init.pdb 

        structure water.pdb 
          number {}
          inside cube -{} -{} -{} {}
        end structure 

        """.format(n_water, box_size/2.0, box_size/2.0, box_size/2.0, box_size)
        , file=fout)

    print('='*80)
    with open('packmol.inp', 'r') as f:
        print(f.read())
    print('='*80)

    result = run_command(opt.packmol )
    with open('packmol.out', 'w') as f:
        print(result, file=f)


    gmxscript.editconf(
        f = "system_init.pdb",
        box = [box_size/10.0, box_size/10.0, box_size/10.0],
        bt = "cubic",
        o = "waterbox.pdb"
    )

    with open('waterbox.top', 'w') as f:
        print(tip4p_fb_itp.format(n_water), file=f)

    with open('minim.mdp', 'w') as f:
        print(minim_mdp.format(half_box, half_box), file=f)
    with open('nvt.mdp', 'w') as f:
        print(nvt_mdp.format(half_box, half_box), file=f)

    gmxscript.grompp(f="minim.mdp", c="waterbox.pdb", o="mini.tpr", p="waterbox.top")
    gmxscript.mdrun(deffnm="mini")
    gmxscript.select(f="mini.trr", s="mini.tpr", on="mini", select="resname SOL")
    gmxscript.grompp(f="nvt.mdp", c="mini.gro", p="waterbox.top", o="nvt.tpr")
    gmxscript.mdrun(deffnm="nvt", v=True)    
    gmxscript.trjconv(f="nvt.trr", s="nvt.tpr", n="mini.ndx", o="final_water_box.g96")
    os.chdir(cwd)

# if not os.path.exists('prep/final_water_box.g96'):
shutil.rmtree('prep', ignore_errors=True)
run_prep_box()



with open('prep/final_water_box.g96', 'r') as f, open(opt.output, 'w')  as fout, open('output.xyz', 'w')  as fxyz:
    print(f"{n_water*3} 0 1", file=fout)
    print(f"{n_water*3}", file=fxyz)
    print("", file=fxyz)
    for line in f:
        if line.startswith("POSITIONRED"):
            for i in range(n_water*4):
                line = next(f)
                if (i+1) % 4 == 0: 
                    continue
                sym = 'H'
                if i%4 == 0:
                    sym = 'O'
                arr = list(map(lambda x:float(x)*10.0, line.split()))
                print ('{:4s} {:15.6f} {:15.6f} {:15.6f}'.format(sym, *arr), file=fout)
                print ('{:4s} {:15.6f} {:15.6f} {:15.6f}'.format(sym, *arr), file=fxyz)
            break
    latt = np.diag([box_size]*3).tolist()
    for i in range(3):
        print ('{:4s} {:15.6f} {:15.6f} {:15.6f}'.format('TV', *latt[i]), file=fout)
    if (opt.subsystem):
        print ("", file=fout)
        dig = len(f"{n_water*3}")
        length = int(80/((dig+2)*3))
        for i in range(n_water):
            for_str = " {{:{}d}} ".format(dig)
            print(for_str.format(i+1)*3, end='', file=fout)
            if (i+1) % length == 0:
                print ("", file=fout)
        print ("", file=fout)
        print ("", file=fout)

# shutil.rmtree('prep', ignore_errors=True)




