#!/usr/bin/env python3
import sys
import pymatgen as mg
import pymatgen.core.surface as sf
import pymatgen.symmetry.analyzer as sa
import numpy as np
import json
import os
import shutil
import utils
import scipy.constants

Hartree2Ev = scipy.constants.physical_constants['Hartree energy in eV'][0]
eV2Joule = scipy.constants.physical_constants['electron volt-joule relationship'][0]

class SurfaceEnergyEvaluator():

    def round_up_to_odd(self, f):
        return int(np.ceil(f) // 2 * 2 + 1)

    def KpointsToDFTB(self, kpoints, filename='in.kpoints'):
        with open(filename, 'w') as f:
            kpts = kpoints.kpts[0]
            print('{} 0 0'.format(kpts[0]), file=f)
            print('0 {} 0'.format(kpts[1]), file=f)
            print('0 0 {}'.format(kpts[2]), file=f)
            if (kpoints.style == kpoints.supported_modes.Gamma):
                print('0 0 0', file=f)
            else:
                sh1 = 0.0
                sh2 = 0.0
                sh3 = 0.0
                if (kpts[0] % 2 == 0):
                    sh1 = 0.5
                if (kpts[1] % 2 == 0):
                    sh2 = 0.5
                if (kpts[2] % 2 == 0):
                    sh3 = 0.5
                print('{} {} {}'.format(sh1, sh2, sh3), file=f)

    def __init__(self):
        self.max_hkl_index = 3
        self.min_slab_thickness = 20
        self.min_vacuum_thickness = 20
    
    def compute_surface_energy(self, root_path, struct, skfile_path, temp_input_path, dftbplus_path):
                
        running_path = os.path.join(root_path, 'surface_energy')
        os.makedirs(running_path, exist_ok=True)

        DFTB_opt_section = """
        !Driver = !LBFGS{
            AppendGeometries = Yes
        }
        """

        sga = sa.SpacegroupAnalyzer(struct)
        conv_bulk = sga.get_conventional_standard_structure()

        print('Generating slabs up to index {}'.format(self.max_hkl_index))

        

        millers = sf.get_symmetrically_distinct_miller_indices(conv_bulk, self.max_hkl_index)
        millers.sort()
        # millers = [(3,2,1)]

        res = []

        for miller in millers:
            gen = sf.SlabGenerator(conv_bulk, miller, self.min_slab_thickness,
                                self.min_vacuum_thickness, center_slab=True, primitive=True,
                                max_normal_search=self.max_hkl_index)
            slabs = gen.get_slabs(symmetrize=False)

            print("%s has %d slabs... " % (miller, len(slabs)))
            slab = slabs[0]

            print('Number of atoms: OUC: {}, SLAB: {}'.format(len(slab.oriented_unit_cell.sites), len(slab.sites)))

            print('Surface area: {} A^2'.format(slab.surface_area))

            hkl_str = '{}{}{}'.format(*slab.miller_index)
            temp_path = os.path.join(running_path, hkl_str)
            os.makedirs(temp_path, exist_ok=True)

            ouc_folder = os.path.join(temp_path, 'ouc')
            slab_folder = os.path.join(temp_path, 'slab')

            os.makedirs(ouc_folder, exist_ok=True)
            os.makedirs(slab_folder, exist_ok=True)

            # filename = os.path.join(slab_folder, 'in.gen')
            # utils.StructureToGen(slab, filename=filename)

            # filename = os.path.join(ouc_folder, 'in.gen')
            # utils.StructureToGen(slab.oriented_unit_cell, filename=filename)

            kptsa = self.round_up_to_odd(50.0/slab.oriented_unit_cell.lattice.a)
            kptsb = self.round_up_to_odd(50.0/slab.oriented_unit_cell.lattice.b)
            kptsc = self.round_up_to_odd(50.0/slab.oriented_unit_cell.lattice.c)
            kpoints_obj_ouc = mg.io.vasp.Kpoints.gamma_automatic(kpts=(kptsa, kptsb, kptsc))
            kpoints_obj_slab = mg.io.vasp.Kpoints.gamma_automatic(kpts=(kptsa, kptsb, 1))

            print('Comuting bulk energy using NOUC ... ', end='')
            l_res = utils.evaluate_dftb(slab.oriented_unit_cell, ouc_folder, skfile_path=skfile_path, temp_input_path=temp_input_path, dftbplus_path=dftbplus_path, kpoints=kpoints_obj_ouc, geom_opt=False)
            ouc_energy = l_res['energy_final']
            print(ouc_energy, ' H')

            print('Comuting surface ... ', end='')
            l_res = utils.evaluate_dftb(slab, slab_folder, skfile_path=skfile_path, temp_input_path=temp_input_path, dftbplus_path=dftbplus_path, kpoints=kpoints_obj_slab, geom_opt=True)
            slab_energy = l_res['energy_final']
            print(slab_energy, ' H')
            surface_energy = (slab_energy - len(slab.sites)*ouc_energy/len(slab.oriented_unit_cell.sites))/(2.0*slab.surface_area)
            surface_energy_Jm2=surface_energy*Hartree2Ev*eV2Joule*10**20
            print('Surface energy: ', surface_energy_Jm2)

            res.append( (miller, surface_energy_Jm2)) 
        return res
            

        # all_hkls = []

        # curdir = os.curdir
        # # os.makedirs('output', exist_ok=True)
        # # os.chdir('output')

        # areas = []

        # for slab in slabs_all:
            
        #     print('Writing slab: {}'.format(hkl_str))

        #     ouc_folder = os.path.join('output', hkl_str, 'ouc')
        #     slab_folder = os.path.join('output', hkl_str, 'slab')
        #     os.makedirs(ouc_folder, exist_ok=True)
        #     os.makedirs(slab_folder, exist_ok=True)

        #     filename = os.path.join(ouc_folder, 'slab_{}_ouc.vasp'.format(hkl_str) )
        #     slab.oriented_unit_cell.to(fmt='poscar', filename=filename)
        #     filename = os.path.join(slab_folder, 'slab_{}.vasp'.format(hkl_str) )
        #     slab.to(fmt='poscar', filename=filename)

        #     filename = os.path.join(ouc_folder, 'in.gen')
        #     StructureToGen(slab.oriented_unit_cell, filename=filename)
        #     filename = os.path.join(slab_folder, 'in.gen')
        #     StructureToGen(slab, filename=filename)

        #     all_hkls.append(hkl_str)
        #     values = {}
        #     values['ouc'] = slab.oriented_unit_cell.as_dict()
        # #    values['slab'] = slab.to_json()
        #     values['area'] = slab.surface_area
        #     values['nat'] = len(slab.sites)
        #     values['shift'] = slab.shift
        #     values['scale_factor'] = slab.scale_factor.tolist()
            
        #     areas.append( (hkl_str, slab.surface_area) )
        
        #     kptsa = round_up_to_odd(50.0/slab.oriented_unit_cell.lattice.a)
        #     kptsb = round_up_to_odd(50.0/slab.oriented_unit_cell.lattice.b)
        #     kptsc = round_up_to_odd(50.0/slab.oriented_unit_cell.lattice.c)
        #     kpoints_obj = mg.io.vasp.Kpoints.gamma_automatic(kpts=(kptsa, kptsb, kptsc))
            
        #     filename = os.path.join(ouc_folder, 'KPOINTS')
        #     kpoints_obj.write_file(filename)
        #     filename = os.path.join(ouc_folder, 'in.kpoints')
        #     KpointsToDFTB(kpoints_obj, filename=filename)


        #     kpoints_obj = mg.io.vasp.Kpoints.gamma_automatic(kpts=(kptsa, kptsb, 1))
        #     filename = os.path.join(slab_folder, 'KPOINTS')
        #     kpoints_obj.write_file(filename)
        #     filename = os.path.join(slab_folder, 'in.kpoints')
        #     KpointsToDFTB(kpoints_obj, filename=filename)
        #     filename = os.path.join(slab_folder, 'area')
        #     with open(filename, 'w') as f:
        #         print(slab.surface_area, file=f)

        #     for temp_file, temp_name in zip(template_inputs_abs, template_inputs):
        #         # print(temp_file, ouc_folder)
        #         shutil.copyfile(temp_file, os.path.join(ouc_folder, temp_name))
        #         shutil.copyfile(temp_file, os.path.join(slab_folder, temp_name))

        #     filename = os.path.join(slab_folder, 'dftb_in.hsd')
        #     with open(filename, 'a') as f:
        #         print(DFTB_opt_section, file=f)
            


        #     metadata[hkl_str] = values

        # metadata['hkls'] = all_hkls

        # # os.chdir(curdir)

        # with open('slab_info.json', 'w') as f:
        #     json.dump(metadata, f, indent=2)

        # areas.sort()
        # for hkl, area in areas:
        #     print('{:5s} {:20.12f}'.format(hkl, area))
