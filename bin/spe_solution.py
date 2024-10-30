#!/usr/bin/env python3

import subprocess
import os
import shutil
import gmxscript
from scipy.constants import N_A
import argparse
import sys
from cmmde_mass import mass
parser = argparse.ArgumentParser(description='Computational Molecular & Material Design Interface')
parser.add_argument('-mt', '--terlarut', type=str, help='Nama file molekul terlarut dalam format .xyz (ditulis tanpa ekstensi).')
parser.add_argument('-ct', '--c_terlarut', type=int, help='Muatan bersih molekul terlarut.',default=0)
parser.add_argument('-conc', '--concentration', type=float, default=1.0, help='Konsentrasi molekul terlarut dalam molal (mol/kg).')
parser.add_argument('-mp', '--pelarut', type=str, help='Nama file molekul pelarut dalam format .xyz (ditulis tanpa ekstensi).')
parser.add_argument('-cp', '--c_pelarut', type=int, help='Muatan bersih molekul pelarut',default=0)
parser.add_argument('-Nump','--NumPelarut', type=int, default=100, help='Jumlah molekul pelarut maksimum dalam sistem larutan. Default = 100.')
parser.add_argument('-cat','--cation',type=str, default='Li',help='Kation yang digunakan dalam sistem elektrolit baterai.')
parser.add_argument('-p', '--pressure', default=1.0, type=float, help='Tekanan dalam satuan bar')
parser.add_argument('-t', '--temperature', default=298.15, type=float, help='Suhu yang digunakan dalam simulasi.')
parser.add_argument('-nequil','--nequil',type=int, default=50000, help='Jumlah step pada saat ekuilibrasi.')
parser.add_argument('-nprod','--nprod',type=int, default=400000, help='Jumlah step pada saat production.')
parser.add_argument('-dt','--dt',type=float, default=0.0005, help='Step simulasi dinamika molekul yang dilakukan.')
parser.add_argument('--packmol', type=str, default='/home/adit/opt/packmol/packmol')
opt=parser.parse_args(sys.argv[1:])

def run_command(command, args=[], input_str=None):
    command_list = [command]
    if (len(args) > 0):
        command_list.extend(args)

    if (input_str is not None):
        res =  subprocess.run(command_list, stdout=subprocess.PIPE, input=input_str.encode(encoding='ascii'))
    else:
        res =  subprocess.run(command_list, stdout=subprocess.PIPE)
    return res.stdout.decode()


# Input file minimisasi energi dalam GROMACS
minim_mdp = """
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 100.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 5         ; Frequency to update the neighbor list and long range forces
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = {}       ; Short-range electrostatic cut-off
rvdw            = {}       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
;define = -DFLEXIBLE
"""
nve_mdp = """
define                  = -DFLEXIBLE
integrator              = md        ; leap-frog integrator
nsteps                  = {}    ; 1 * 500000 = 50 ps
dt                      = 0.001     ; 0.5 fs
nstxout                 = 100       ; save coordinates every 1.0 ps
nstenergy               = 100       ; save energies every 1.0 ps
nstlog                  = 100       ; update log file every 1.0 ps

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
; Temperature coupling is off
tcoupl                  = no             ; no temperature coupling in NVE
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVE
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no       ; Velocity generation is off
"""

# Input file simulasi menggunakan ensembel kanonik (NVT) dalam GROMACS
nvt_mdp="""
define                  = -DFLEXIBLE
integrator              = md        ; leap-frog integrator
nsteps                  = {}    ; 1 * 500000 = 50 ps
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
ref_t                   = {}              ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = {}       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
"""

# Input file simulasi isoterm-isobarik (NPT) dalam GROMACS. Harapannya, setelah dilakukan NPT, akan diperoleh sistem dengan massa jenis yang optimal.
npt_mdp = """
title		= NPT equilibration 
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= {}		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = {}		; short-range electrostatic cutoff (in nm)
rvdw		    = {}		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps     = System   ; two coupling groups - more accurate
tau_t		= 0.1	  	        ; time constant, in ps
ref_t		= {} 	  	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Berendsen	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = {}		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off 
"""

nvep_mdp = """
title       = NVT Production
; Run parameters
integrator  = md        ; leap-frog integrator
nsteps      = {}        ; 
dt          = {}     ; in ps
; Output control
nstxout     = 500       ; save coordinates every 1.0 ps
nstvout     = 500       ; save velocities every 1.0 ps
nstenergy   = 500       ; save energies every 1.0 ps
nstlog      = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes        ; first dynamics run
;constraint_algorithm    = lincs     ; holonomic constraints
;constraints             = h-bonds   ; bonds involving H are constrained
;lincs_iter              = 1         ; accuracy of LINCS
;lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type         = grid      ; search neighboring grid cells
nstlist         = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb        = {}        ; short-range electrostatic cutoff (in nm)
rvdw            = {}        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4         ; cubic interpolation
fourierspacing  = 0.16      ; grid spacing for FFT  
; Temperature coupling is off
;tcoupl = no
tcoupl      = V-rescale             ; Tcoupling off in NVE
tc-grps     = System   ; two coupling groups - more accurate
tau_t   = 0.1             ; time constant, in ps
ref_t   = {}              ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl            = Berendsen     ; Pressure coupling on in NPT
pcoupltype          = isotropic             ; uniform scaling of box vectors
tau_p           = 2.0               ; time constant, in ps
ref_p           = {}                ; reference pressure, in bar
compressibility     = 4.5e-5              ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel     = no        ; Velocity generation is off
"""
def make_solution(terlarut,c_terlarut,pelarut,c_pelarut,N_pelarut,concentration,cation, temperature,pressure,packmol):

  if temperature:
    print (f'Suhu: {temperature:12.4f} K')

# Penyiapan box molekul
  os.chdir('{}'.format(terlarut))
# Memanggil antechamber untuk membuatkan file mol2 yang mengandung informasi muatan bersih atom
  os.system("obabel cmmd.xyz -O {}.pdb".format(terlarut))
  with open('{}.pdb'.format(terlarut),'r') as f:
    fdata = f.read()
  fdata = fdata.replace('UNL','SOL')
  fdata = fdata.replace('CONECT', '#CONECT')
  with open('{}.pdb'.format(terlarut),'w') as f:
    f.write(fdata) 
  os.system("antechamber -j 5 -at gaff -dr no -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 -c bcc -s 2 -nc {}".format(terlarut,terlarut,c_terlarut))
# Mengidentifikasi parameter yang tidak terdefinisikan dalam GAFF
  os.system("parmchk2 -i {}.mol2 -f mol2 -o {}.frcmod".format(terlarut,terlarut))
  with open("tleap.in", 'w') as fout:
    print("""
    source leaprc.gaff
        SOL = loadmol2 {}.mol2
        check SOL
        loadamberparams {}.frcmod
        saveoff SOL {}.lib
        saveamberparm SOL {}.parmtop {}.inpcrd
        savepdb SOL {}_gaff.pdb
        quit
    """.format(terlarut,terlarut,terlarut,terlarut,terlarut,terlarut), file=fout)
  os.system("tleap -f tleap.in")

  os.chdir('../{}'.format(pelarut))

# Penyiapan molekul pelarut menggunakan GAFF
# Memanggil antechamber untuk membuatkan file mol2 yang mengandung informasi muatan bersih atom
  os.system("obabel cmmd.xyz -O {}.pdb".format(pelarut))
  with open('{}.pdb'.format(pelarut),'r') as f:
    fdata = f.read()
  fdata = fdata.replace('UNL','SLV')
  fdata = fdata.replace('CONECT', '#CONECT')
  with open('{}.pdb'.format(pelarut),'w') as f:
    f.write(fdata)
  os.system("antechamber -j 5 -at gaff -dr no -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 -c bcc -s 2 -nc {}".format(pelarut,pelarut,c_pelarut))
# Mengidentifikasi parameter yang tidak terdefinisikan dalam GAFF
  os.system("parmchk2 -i {}.mol2 -f mol2 -o {}.frcmod".format(pelarut,pelarut))

  with open("tleap.in", 'w') as fout:
    print("""
    source leaprc.gaff
        SLV = loadmol2 {}.mol2
        check SLV
        loadamberparams {}.frcmod
        saveoff SLV {}.lib
        saveamberparm SLV {}.parmtop {}.inpcrd
        savepdb SLV {}_gaff.pdb
        quit
  """.format(pelarut,pelarut,pelarut,pelarut,pelarut,pelarut), file=fout)
  os.system("tleap -f tleap.in")
  os.system("mkdir ../{}".format(cation))
  os.chdir('../{}'.format(cation))
  with open('{}.mol2'.format(cation),'w') as f:
    print("""@<TRIPOS>MOLECULE
ION
 1 0 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 {}+          0.0000    0.0000    0.0000 {}+      1  ion        1.0000
@<TRIPOS>BOND
    """.format(cation,cation),file=f)
  os.chdir('../')
# Penyiapan box
  Mr2 = 0
  elements = []
  with open ("{}/cmmd.xyz".format(pelarut), 'r') as f:
    next(f)
    next(f)
    for line in f:
      arr = line.split()
      elements.append(arr[0])

  for ele in elements:
    Mr2 += float(mass(ele))
    
  m_pelarut = N_pelarut*(Mr2/N_A)/1000 # dalam kg


  mol_ions = concentration*m_pelarut
  N_ions = round(mol_ions*N_A)
  V_box = m_pelarut*1000 # Dalam satuan cm3
# Panjang rusuk kubus dalam satuan cm:
  box_size = V_box**(1/3) 
# Panjang rusuk kubus dalam satuan Angstroem:
  cmToAng = 1e8
  box_size = box_size * cmToAng + 5
# Hentikan perhitungan jika ukuran kubus terlalu kecil
  half_box = min(1.0, box_size/20.0-0.1)
  print(half_box)

  def run_prep_box():
      cwd = os.path.abspath(os.curdir)
      print('Membuat sistem larutan menggunakan PACKMOL')
      os.makedirs('SistemLarutan', exist_ok=True)
      os.chdir('SistemLarutan')


      with open('packmol.inp', 'w') as fout:
          print("""tolerance 2.5
filetype pdb
output system_init.pdb 

structure ../{}/{}_gaff.pdb 
    number {}
    inside cube -{} -{} -{} {}
    resnumbers 3
  end structure 

structure ../{}/{}_gaff.pdb 
    number {}
    inside cube -{} -{} -{} {}
    resnumbers 3
  end structure """.format(pelarut, pelarut, N_pelarut, box_size/2.0, box_size/2.0, box_size/2.0, box_size, terlarut, terlarut, N_ions, box_size/2.0, box_size/2.0, box_size/2.0, box_size), file=fout)
      print('='*80)

      with open ('tleap.in', 'w') as fout:
      	print("""source leaprc.gaff
loadamberparams /home/adit/miniconda3/dat/leap/parm/frcmod.ions1lm_iod
loadoff ../{}/{}.lib
loadoff ../{}/{}.lib
loadamberparams ../{}/{}.frcmod
loadamberparams ../{}/{}.frcmod
ION = loadmol2 ../{}/{}.mol2
SYSTEM = loadpdb system_init.pdb
addions SYSTEM ION 0
list
saveamberparm SYSTEM system.prmtop system.inpcrd
savepdb SYSTEM system.pdb
quit""".format(terlarut,terlarut,pelarut,pelarut,terlarut,terlarut,pelarut,pelarut,cation,cation), file=fout)

      with open('packmol.inp', 'r') as f:
          print(f.read())
      print('='*80)

      result = run_command(packmol)
      with open('packmol.out', 'w') as f:
          print(result, file=f)

      os.system("tleap -f tleap.in")
      os.system("acpype -p system.prmtop -x system.inpcrd")

      gmxscript.editconf(
          f = "system.pdb",
          box = [box_size/10.0, box_size/10.0, box_size/10.0],
          bt = "cubic",
          o = "system_box.pdb"
      )

  # Menuliskan semua input file yang dibutuhkan, meliputi simulasi NVE, NVT, dan NPT
      nequil = opt.nequil
      nprod = opt.nprod

      with open('minim.mdp', 'w') as f:
          print(minim_mdp.format(half_box, half_box), file=f)
      with open('nve.mdp', 'w') as f:
          print(nve_mdp.format(nequil,half_box,half_box), file=f)
      with open('nvt.mdp', 'w') as f:
          print(nvt_mdp.format(nequil,half_box, half_box,temperature,temperature), file=f)
      with open('npt.mdp', 'w') as f:
        print(npt_mdp.format(nequil,half_box, half_box,temperature,pressure), file=f)
      with open('nvep.mdp', 'w') as f:
          print(nvep_mdp.format(nprod,opt.dt,half_box,half_box,temperature,pressure), file=f)
      gmxscript.grompp(f="minim.mdp", c="system_box.pdb", o="mini.tpr", p="system_GMX.top", maxwarn=2)
      gmxscript.mdrun(deffnm="mini",ntmpi=1)
      gmxscript.select(f="mini.trr", s="mini.tpr", on="mini", select="System")
      gmxscript.grompp(f="nve.mdp", c="mini.gro", p="system_GMX.top", o="nve.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="nve", v=True, ntmpi=1)
      gmxscript.grompp(f="nvt.mdp", c="nve.gro", p="system_GMX.top", o="nvt.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="nvt", v=True, ntmpi=1)   
      gmxscript.grompp(f="npt.mdp", c="nvt.gro", p="system_GMX.top", t="nvt.cpt",o="npt.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="npt", v=True, ntmpi=1)
      gmxscript.grompp(f="nvep.mdp", c="npt.gro", p="system_GMX.top", o="nvep.tpr",maxwarn=2)
      gmxscript.mdrun(deffnm="nvep", v=True, ntmpi=1)   
      gmxscript.trjconv(f="nvep.trr", s="nvep.tpr", n="mini.ndx", o="finalsystem.pdb")
      os.system("obabel finalsystem.pdb -O finalsystem.xyz")
      with open('finalsystem.xyz', 'r') as f:
        fdata = f.read()
      fdata = fdata.replace('*','{}'.format(opt.cation))
      with open('finalsystem.pdb', 'r') as f:
        tv1 = 0
        tv2 = 0
        tv3 = 0
        for line in f:
          if 'CRYST' in line:
            arr = line.split()
            tv1 += float(arr[1])
            tv2 += float(arr[2])
            tv3 += float(arr[3])
            break
      with open('finalsystem.xyz'.format(terlarut),'w') as fout:
        fout.write(fdata)
        print("TV {} 0.00 0.00".format(tv1),file=fout)
        print("TV 0.00 {} 0.00".format(tv2), file=fout)
        print("TV 0.00 0.00 {}".format(tv3), file=fout)
      os.system("/home/adit/opt/gromacs-2021.3/build/bin/gmx msd -f nvep.trr -s nvep.tpr -o msd_SOL.xvg")
      os.system("/home/adit/opt/gromacs-2021.3/build/bin/gmx msd -f nvep.trr -s nvep.tpr -o msd_ION.xvg")
      os.system("/home/adit/opt/gromacs-2021.3/build/bin/gmx msd -f nvep.trr -s nvep.tpr -o msd_SLV.xvg")
      os.chdir(cwd)
    
  shutil.rmtree('SistemLarutan', ignore_errors=True)
  run_prep_box()
# Eksekusi simulasi
make_solution(opt.terlarut,opt.c_terlarut,opt.pelarut,opt.c_pelarut,opt.NumPelarut,opt.concentration,opt.cation,opt.temperature,opt.pressure,opt.packmol)




