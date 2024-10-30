#!/usr/bin/env python3

import sys
import new_rep
import genformat
import spline_rep
import pymatgen
import math
import numpy as np
import scipy.constants


def evaluate_crystal(reps: dict, struct: pymatgen.core.IStructure):
    force_factor = 51.42208619083232
    max_cutoffs = {}
    for k, v in reps.items():
        elems = k.split('-')
        max_cutoffs[elems[0]] = max(max_cutoffs.get(elems[0],0.0), v.cutoff)
        max_cutoffs[elems[1]] = max(max_cutoffs.get(elems[1],0.0), v.cutoff)
    energy = 0.0
    forces = np.zeros((len(struct.sites), 3))
    for ind, site in enumerate(struct.sites):
        neighbors = struct.get_neighbors(site, max_cutoffs[str(site.specie)]/pymatgen.core.units.ang_to_bohr, include_index=True)
        for n_site, n_r, n_index in neighbors:   
            rep_name = '{}-{}'.format(str(site.specie), str(n_site.specie))
            r_in_bohr = n_r*pymatgen.core.units.ang_to_bohr           
            derv0 = reps[rep_name].eval(r_in_bohr)

            if ind == n_index:
                energy += 0.5*derv0
            else:
                if (ind < n_index):    
                    energy += derv0  
                derv1 = reps[rep_name].eval(r_in_bohr, 1)
                vec = n_site.coords - site.coords    
                forces[ind] += (-derv1) * vec / n_r
                # forces[n_index] -= (-derv1) * vec / n_r
                
    return energy*scipy.constants.physical_constants['atomic unit of electric potential'][0], forces*force_factor
    # return energy, forces

def evaluate_molecule(reps: dict, struct: pymatgen.core.IMolecule):
    max_cutoffs = {}
    for k, v in reps.items():
        elems = k.split('-')
        max_cutoffs[elems[0]] = max(max_cutoffs.get(elems[0],0.0), v.cutoff)
        max_cutoffs[elems[1]] = max(max_cutoffs.get(elems[1],0.0), v.cutoff)
    energy = 0.0
    forces = np.zeros((len(struct.sites), 3))
    for ind, site in enumerate(struct.sites):
        for ind2, site2 in enumerate(struct.sites):
            if (ind<=ind2):
                continue
            vec = site2.coords - site.coords
            n_r = math.sqrt(sum(vec*vec))
            rep_name = '{}-{}'.format(str(site.specie), str(site2.specie))
            r_in_bohr = n_r*pymatgen.core.units.ang_to_bohr

            derv0 = reps[rep_name].eval(r_in_bohr)
            derv1 = reps[rep_name].eval(r_in_bohr, 1)

            energy += derv0
            forces[ind] += (-derv1) * vec / n_r
            forces[ind2] -= (-derv1) * vec / n_r

    return energy, forces

def evaluate(reps, struct):
    if (isinstance(struct, pymatgen.core.IStructure)):
        return evaluate_crystal(reps, struct)
    elif (isinstance(struct, pymatgen.core.IMolecule)):
        return evaluate_molecule(reps, struct)


if __name__ == "__main__":  
    # coeffs = [ 1.33354733e-01, -2.75117922e-01,  1.78036236e-03,  3.98357310e-01,
    # -1.28598691e-01, -6.00153847e-01,  9.10383626e-01, -6.48505490e-01,
    # 2.77216939e-01, -7.49863115e-02,  1.26340420e-02, -1.21510402e-03,
    # 5.11081519e-05]
    
    # rep = new_rep.PolynomialRepulsive(6.5, coeffs)

    spline_rep = spline_rep.RepulsivePotenial.from_file('Pd-Pd.skf')

    reps = {}
    reps['Pd-Pd'] = spline_rep

    input_geom_file = sys.argv[1]

    geom = genformat.loadgen(input_geom_file)

    energy, forces = evaluate(reps, geom)
    
    print('{:20.12F}'.format(energy))
    for force in forces:
        print('{:20.12F} {:20.12F} {:20.12F}'.format(*force))



