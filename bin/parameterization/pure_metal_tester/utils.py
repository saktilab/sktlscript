#!/usr/bin/env python3

import shutil
import pymatgen
import os
import subprocess

def StructureToGen(struct, filename='in.gen'):
    nat = len(struct.sites)
    symbols = [str(x.specie) for x in struct.sites]
    uniq_symbols = list(set(symbols))
    with open(filename, 'w') as f:
        print('{} F'.format(nat), file=f)
        print(' '.join(uniq_symbols),file=f)
        for ind, sym in enumerate(symbols):
            print ('{:6d} {:3d} {:20.12f} {:20.12f} {:20.12f}'.format(
                ind+1, uniq_symbols.index(sym)+1, *struct.sites[ind].frac_coords
            ),file=f)
        print('0.0 0.0 0.0', file=f)
        for i in range(3):
            print('{:20.12f} {:20.12f} {:20.12f}'.format(*struct.lattice.matrix[i]), file=f)

def load_gen(filename):
    with open(filename, 'r') as f:
        line = next(f)
        arr = line.split()
        gentype = arr[1]
        nat = int(arr[0])
        elements = next(f).split()
        species = []
        coords = []
        for i in range(nat):
            line = next(f)
            arr = line.split()
            species.append(elements[int(arr[1])-1])
            coords.append(list(map(float, arr[2:5])))
        if gentype == 'C':
            return pymatgen.core.Molecule(species=species, coords=coords)
        else:
            line = next(f)
            tv = []
            for i in range(3):
                line = next(f)
                tv.append(list(map(float, line.split())))

            cart_coords = True
            if gentype == 'F':
                cart_coords = False
            return pymatgen.core.Structure(pymatgen.Lattice(tv), species=species, coords=coords, coords_are_cartesian=cart_coords)




class DummyBlock:
    def __init__(self, path):
        self.path = path
    def __enter__(self):
        path = self.path
        os.makedirs(path, exist_ok=True)
        return path
        
    def __exit__(self, exc_type, exc_value, traceback):
        pass

def evaluate_dftb(struct, path, skfile_path, temp_input_path, dftbplus_path, kpoints=None, geom_opt=False, opt_lattice=False, opt_fixangles=False):

    if (kpoints == None):
        kpoints = pymatgen.io.vasp.Kpoints.automatic_gamma_density(struct, 8000)

    l_curdir = os.path.abspath(os.curdir)
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    
    # with tempfile.TemporaryDirectory(prefix='/ram/dftb/') as tmpdir:
    with DummyBlock(path) as tmpdir:
        os.chdir(tmpdir)
        # shutil.copyfile(os.path.join(temp_input_path, 'dftbinp', 'crystal_energy.hsd'), 'dftb_in.hsd')
        shutil.copyfile(temp_input_path, 'dftb_in.hsd')
        with open('dftb_in.hsd', 'a') as f:
            print("""
Geometry = GenFormat{
<<< in.gen
}
""",file=f)
            sa = 0.0
            sb = 0.0
            sc = 0.0
            if (kpoints.style == pymatgen.io.vasp.Kpoints.supported_modes.Monkhorst):           
                if (kpoints.kpts[0][0] % 2):
                    sa = 0.5
                
                if (kpoints.kpts[0][1] % 2):
                    sb = 0.5
                
                if (kpoints.kpts[0][2] % 2):
                    sc = 0.5
            print("""
+Hamiltonian = +DFTB{{
    !SlaterKosterFiles = !Type2FileNames{{
        prefix = '{}/'
        suffix = '.skf'
        separator = '-'
    }}
    !KPointsAndWeights = SupercellFolding{{
        {} 0 0
        0 {} 0
        0 0 {}
        {} {} {}
    }}
}}
""".format(skfile_path, *kpoints.kpts[0], sa, sb, sc), file=f)
            if (geom_opt):
                opt_strs = []
                if (opt_lattice):
                    opt_strs.append('    LatticeOpt = Yes')
                if (opt_fixangles):
                    opt_strs.append('    FixAngles = Yes')
                print("""
    !Driver = !LBFGS{{
        {}
    }}
    """.format('\n'.join(opt_strs)),file=f)

            else:
                print("""
    !Driver = {{}}
    """,file=f)


        StructureToGen(struct)
        res = subprocess.run([dftbplus_path], stdout=subprocess.PIPE)
        # print(res)
        with open('dftb.out', 'w') as fout:
            print(res.stdout.decode("utf-8"), file=fout)        

        try:
            with open('detailed.out', 'r') as f:
                for line in f:
                    if 'Total Mermin free energy:' in line:
                        energy = float(line.split()[4])
        except Exception as ex:
            print(os.path.abspath(os.curdir))
            raise ex

        res = {}
        res['energy_final'] = energy

        if (geom_opt):
            optimized_geom = load_gen('geo_end.gen')
            res['geom_final'] = optimized_geom
        
        
        

        os.chdir(l_curdir)
        return res

def evaluate_dftb_input(path, skfile_path, temp_input_path, dftbplus_path):

    l_curdir = os.path.abspath(os.curdir)
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    
    # with tempfile.TemporaryDirectory(prefix='/ram/dftb/') as tmpdir:
    with DummyBlock(path) as tmpdir:
        os.chdir(tmpdir)
        # shutil.copyfile(os.path.join(temp_input_path, 'dftbinp', 'crystal_energy.hsd'), 'dftb_in.hsd')
        shutil.copyfile(temp_input_path, 'dftb_in.hsd')
        with open('dftb_in.hsd', 'a') as f:
            print("""
+Hamiltonian = +DFTB{{
    !SlaterKosterFiles = !Type2FileNames{{
        prefix = '{}/'
        suffix = '.skf'
        separator = '-'
    }}
}}
""".format(skfile_path), file=f)
            
        res = subprocess.run([dftbplus_path], stdout=subprocess.PIPE)
        # print(res)
        with open('dftb.out', 'w') as fout:
            print(res.stdout.decode("utf-8"), file=fout)        

        try:
            with open('detailed.out', 'r') as f:
                for line in f:
                    if 'Total Mermin free energy:' in line:
                        energy = float(line.split()[4])
        except Exception as ex:
            print(os.path.abspath(os.curdir))
            raise ex

        res = {}
        res['energy_final'] = energy    
        
        os.chdir(l_curdir)
        return res