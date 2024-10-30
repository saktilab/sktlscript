#!/usr/bin/env python3

import sys
import os
import argparse
import gmxscript

parser = argparse.ArgumentParser(description='Ligand parameter generator')
parser.add_argument('-l', '--ligand', type=str, help='ligand structure in pdb format')
parser.add_argument('-c', '--charge', type=int, help='Net charge of the ligand', default=0)
parser.add_argument('-p', '--protein', type=str, help='protein structure in pdb format')
opt = parser.parse_args(sys.argv[1:])


####GROMACS INPUT FILES####
# Dummy input for ion generation
with open('ions.mdp', 'w') as f:
    print("""
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		; Cut-off for making neighbor list (short range forces)
coulombtype	    = cutoff	; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		; long range electrostatic cut-off
rvdw		    = 1.2		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
           """, file=f)

# Energy minimization
em_mdp = """title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
"""

# Molecular dynamics
md_mdp = """integrator  = md
nsteps      = 0
dt          = 0.002

nstxtcout   = 1000
nstxout     = 1000

constraints = all-bonds
pbc         = xyz

ns_type     = grid
rlist       = 1.4

coulombtype = PME
rcoulomb    = 1.4

vdwtype     = cut-off
rvdw        = 1.4

tcoupl      = V-rescale
tc_grps     = Protein Non-protein
tau_t       = 0.1 0.1
ref_t       = 0 0
"""

ligand = opt.ligand
charge = opt.charge
protein = opt.protein

os.system("antechamber -i {} -fi pdb -o drug.mol2 -fo mol2 -c gas -s 2 -nc {}".format(ligand,charge))

os.system("parmchk2 -i drug.mol2 -f mol2 -o drug.frcmod")

with open("tleap.in", 'w') as fout:
  print("""
  source leaprc.gaff
      LIG = loadmol2 drug.mol2
      check LIG
      loadamberparams drug.frcmod
      saveoff LIG drug.lib
      saveamberparm LIG drug.prmtop drug.inpcrd
      savepdb LIG drug_gaff.pdb
      quit
  """, file=fout)

os.system("tleap -f tleap.in")

os.system("acpype -p drug.prmtop -x drug.inpcrd")

with open("LIG_GMX.top", 'r') as f:
    lines = f.readlines()
    lines = lines[:-3]
with open("LIG_GMX.top", 'w') as f:
    for i, line in enumerate(lines):
        if i not in [2,3,4]:
            f.write(line)

gmxscript.pdb2gmx(
        f = "{}".format(protein),
        ff = "amber99sb",
        water = "spc",
        o = "protein.gro",
        p = "topol.top",
        ignh = "yes"
        )

with open ("protein.gro", 'r') as f:
    next(f)
    nPro = int(next(f))
    pro_atoms = list(range(0,nPro))
    comp = [ ]
    for i, line in enumerate(f):
        if i in pro_atoms:
            comp.append(line.strip("\n"))
        elif i > nPro:
            break
    for line in f:
        pass
    cell = line.strip("\n")

gmxscript.editconf(
        f = "drug_gaff.pdb",
        o = "drug_gaff.gro"
        )


with open("drug_gaff.gro", 'r') as f:
    next(f)
    nLig = int(next(f))
    lig_atoms = list(range(0,nLig))
    for i, line in enumerate(f):
        if i in lig_atoms:
            comp.append(line.strip("\n"))
        elif i > nLig:
            break

with open("complex.gro", 'w') as fout:
    print("Complex.gro", file=fout)
    print("{}".format(nPro+nLig), file=fout)
    for i in comp:
        print(i, file=fout)
    print(cell, file=fout)


with open("topol.top", 'r+') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('[ moleculetype ]'):
            lines[i-1] = lines[i-1].strip() + '\n; Include ligand topology\n' + '#include "LIG_GMX.top"\n\n'
        if line.startswith('[ molecules ]'):
            lines[i+2] = lines[i+2].strip() + '\nLIG                 1\n'
    f.seek(0)
    for line in lines:
        f.write(line)

gmxscript.editconf(
        f = "complex.gro",
        bt = "dodecahedron",
        d = 1.4,
        o = "complex_box.gro",
        c = "yes"
        )


gmxscript.solvate(
        cp = "complex_box.gro",
        cs = "spc216.gro",
        p = "topol.top",
        o = "solv.gro"
        )

gmxscript.grompp(
        f = 'ions.mdp',
        c = 'solv.gro',
        o = 'ions.tpr',
        p = 'topol.top'
        )

gmxscript.genion(
        s = 'ions.tpr',
        o = 'solv_ions.gro',
        p = 'topol.top',
        pname = 'NA',
        nname = 'CL',
        neutral = 'yes',
        stdin = """
        SOL
        """
        )

with open("em.mdp", 'w') as fout:
        print(em_mdp, file=fout)
# steepest descent energy minimization
gmxscript.grompp(
        f = 'em.mdp',
        c = 'solv_ions.gro',
        p = 'topol.top',
        o = 'em.tpr'
        )

gmxscript.mdrun(deffnm='em')

 # Conjugate gradient energy minimization
gmxscript.grompp(
        f = MDP['cg.mdp', {
            'integrator': 'cg',
            'emtol'     : 1.0,
            'nsteps'    : 500,
        }],
        c = 'sd.gro',
        o = 'cg.tpr',
        p = 'topol.top'
    )
gmxscript.mdrun(deffnm = 'cg')

# Position restrained molecular dynamics
# gmxscript.grompp(
#         f = MDP['md.mdp', {
#             'steps'       : 250000,
#             'define'      : '-DPOSRES',
#             'continuation': False,
#             'gen_vel'     : True
#         }],
#         c = 'cg.gro',
#         o = 'mdposres.tpr',
#         p = 'topol.top'
#     )
# gmxscript.mdrun(deffnm = 'mdposres')
