#!/usr/bin/env python3
import pymatgen as mg
import pymatgen.core.units
import math
import numpy as np
import io
import textwrap

from sklearn import linear_model

import utils

import rep_spline
import rep_basis

class RepulsiveModelCollection:
    """ General class for the collection of the repulsive potential model

    """

    def __init__(self, rep_names, rep_models):
        self.rep_names = []
        self.rep_models = []
        self.resulting_repulsives = []
        self.all_names = {}
        self.shifts = {}
        self.max_cutoff = 0.0

        s_rep_names = []
        s_rep_model = []
        for name, model in zip(rep_names, rep_models):
            elems = name.split('-')
            sorted_name = '{}-{}'.format(*sorted(elems))
            s_rep_names.append(sorted_name)
            s_rep_model.append(model)

        shift = 0
        for name, model in sorted(zip(s_rep_names, s_rep_model), key=lambda x: x[0]):
            if (name not in self.all_names):
                self.all_names[name] = model
                elems = name.split('-')
                reversed_name = '{}-{}'.format(*reversed(elems))
                self.all_names[reversed_name] = model
                self.rep_names.append(name)
                self.rep_models.append(model)
                self.shifts[name] = shift
                self.shifts[reversed_name] = shift
                self.max_cutoff = max(self.max_cutoff, model.cutoff)
                shift += model.num_variables
        
        self.total_num_variables = sum([x.num_variables for x  in self.rep_models])

    def description(self):
        fout = io.StringIO()
        print('Potential type: {}'.format(self.rep_models[0].description()), file=fout)
        print('Number of repulsive potentials: {}'.format(len(self.rep_models)), file=fout)
        print('Number of variables: {}'.format(self.get_total_num_variables()), file=fout)
        
        print('Repulsive potenial description:', file=fout)
        for rep_name, rep_model in zip(self.rep_names, self.rep_models):
            print('  {}'.format(rep_name), file=fout)
            print(textwrap.indent(str(rep_model), "    "), file=fout)

        print("", file=fout)
        return fout.getvalue()

    def get_shift(self, name):       
        return self.shifts.get(name, -1)

    def get_num_variables(self, name):
        return self.all_names[name].num_variables        
    
    def get_total_num_variables(self):
        return self.total_num_variables

    def get_row(self, name, r, derv=0):
        return self.all_names[name].get_row(r, derv)
    
class EquationDescription:
    def __init__(self, equ_type: str, equation_item, equ_block, ref_block):
        self.equ_type = equ_type
        self.equation_item = equation_item
        self.equ_block = equ_block
        self.ref_block = ref_block

class RepulsivePotentialEquationBuilder():
    def __init__(self, evaluated_erepfit_input):

        self.evaluated_erepfit_input = evaluated_erepfit_input

        #building a haspmap for easy reference.
        self.uuidmap = {}
        for system in self.evaluated_erepfit_input.systems:
            self.uuidmap[system.uuid] = system

        #equations for each of systems        
        self.equation_descriptions = []

        #load external repulsive potentials
        self.external_repulsive_potentials = { k: rep_spline.RepulsivePotenial.from_file(v) for k, v in self.evaluated_erepfit_input.external_repulsive_potentials.items()}
    
    def get_dist_matrix(self, system, cutoff):
        """ Get the distance and displacement matrix in bohr for all pairs
        """ 
        res = {}
        mole = system.geometry.get_coordinate_struct()
        cutoff_ang = cutoff/mg.core.units.ang_to_bohr
        if (isinstance(mole, mg.Structure)):
            for ind, site in enumerate(mole.sites):
                neighbors = mole.get_neighbors(site, cutoff_ang, include_index=True)
                for n_site, n_r, n_index in neighbors:
                    r_in_bohr = n_r*mg.core.units.ang_to_bohr
                    vec = ((site.coords - n_site.coords)/n_r)
                    if ( (ind, n_index) in res):
                        res[(ind, n_index)].append( (r_in_bohr, vec) )
                    else:
                        res[(ind, n_index)] = [ (r_in_bohr, vec) ]
        else:
            for ind, site in enumerate(mole.sites):
                for ind2, site2 in enumerate(mole.sites):
                    if (ind == ind2):
                        continue
                    vec = site.coords - site2.coords
                    n_r = math.sqrt(sum(vec*vec))
                    r_in_bohr = n_r*mg.core.units.ang_to_bohr
                    p_vec = vec/n_r
                    if (r_in_bohr > cutoff):
                        continue
                    if ( (ind, ind2) in res):
                        res[(ind, ind2)].append( (r_in_bohr, p_vec) )
                    else:
                        res[(ind, ind2)] = [ (r_in_bohr, p_vec) ]
        return res

    def make_energy_block(self, system, dist_matrix, modelcollection, external_reps):
        num_rows = 1
        equ_block = np.zeros( (num_rows, modelcollection.get_total_num_variables()) )
        ref_block = np.zeros( num_rows )
        
        for k, dis in dist_matrix.items():
            ind1, ind2 = k
            rep_name = '{}-{}'.format(system.geometry.coordinates[ind1][0], system.geometry.coordinates[ind2][0])
            shift = modelcollection.get_shift(rep_name)
            if (shift > -1):
                for r, _ in dis:
                    row = modelcollection.get_row(rep_name, r, 0)*0.5
                    lens = len(row)
                    equ_block[0, shift:shift+lens] += row
            elif (rep_name in external_reps):
                for r, _ in dis:
                    ene = external_reps[rep_name].eval(r)
                # print('Add external: {} {} {}'.format(rep_name, r_in_bohr, ene))
                    ref_block -= ene
            else:
                raise RuntimeError('Rep not found: {}'.format(rep_name))
        return equ_block, ref_block
    
    def make_force_block(self, system, dist_matrix, modelcollection, external_reps, filtered_atoms=None):
        
        construct_indexes = set(list(range(len(system.geometry.coordinates))))
        if filtered_atoms is not None:
            construct_indexes = set([ind for ind, x in enumerate(system.geometry.coordinates) if x[0] in filtered_atoms])

        num_rows = 3*len(construct_indexes)
        if num_rows == 0:
            return [], []

        equ_block = np.zeros( (num_rows, modelcollection.get_total_num_variables()))
        ref_block = np.zeros( num_rows )
        
        construct_indexes_index = list(construct_indexes)
        for k, dis in dist_matrix.items():
            ind1, ind2 = k
            if (ind1 not in construct_indexes):
                continue
            ind = construct_indexes_index.index(ind1)
            rep_name = '{}-{}'.format(system.geometry.coordinates[ind1][0], system.geometry.coordinates[ind2][0])
            shift = modelcollection.get_shift(rep_name)
            if (shift > -1):
                for r, vec in dis:
                    row = modelcollection.get_row(rep_name, r, 1)
                    lens = len(row)
                    for xyz in range(3):
                        equ_block[ ind*3+xyz , shift:shift+lens ] += (-row*vec[xyz])
                    
            elif (rep_name in external_reps):
                for r, vec in dis:
                    values = external_reps[rep_name].eval(r, 1)
                    for xyz in range(3):
                        ref_block[ind*3+xyz] -= -values*ref_block[xyz]
            else:
                raise RuntimeError('Rep not found')
            
        return equ_block, construct_indexes_index, ref_block

    def build_linear_equations(self, modelcollection: RepulsiveModelCollection):
        
        self.modelcollection = modelcollection

        # get the maximal cutoff of all repulsive potentials
        max_cutoff = max([x.cutoff for x in self.modelcollection.rep_models])
        if (len(self.external_repulsive_potentials) > 0):
            max_cutoff = max(max_cutoff, max([v.cutoff for k, v in self.external_repulsive_potentials.items()]))

        dist_matrix_uuids = {}
        
        #atomization energy equations
        self._build_energy_equations(dist_matrix_uuids, max_cutoff)

        #force equations
        self._build_force_equations(dist_matrix_uuids, max_cutoff)
    
        #additional equations
        self._build_additional_equations()
        
        #reaction equations
        self._build_reaction_equations(dist_matrix_uuids, max_cutoff)

    def _build_reaction_equations(self, dist_matrix_uuids, max_cutoff):
        for reaction_item in self.evaluated_erepfit_input.equations['reaction']:
            react_energy = utils.get_energy_in_au(reaction_item)

            if (reaction_item.weight < 1.0e-4):
                continue

            final_equ_block = np.zeros( (1, self.modelcollection.get_total_num_variables()) )
            final_ref_block = np.zeros( 1 )

            elec_reaction_energy = 0.0
            for item in reaction_item.products:
                system = self.uuidmap[item.uuid]
                if item.uuid in dist_matrix_uuids:
                    dist_matrix = dist_matrix_uuids[item.uuid]
                else:
                    dist_matrix = self.get_dist_matrix(system, max_cutoff)
                    dist_matrix_uuids[item.uuid] = dist_matrix

                equ_block, ref_block = self.make_energy_block(system, dist_matrix, self.modelcollection, self.external_repulsive_potentials)
                final_equ_block += equ_block*item.coefficient
                final_ref_block += ref_block*item.coefficient
                elec_reaction_energy += system.elec_data['elec_energy']*item.coefficient

            for item in reaction_item.reactants:
                system = self.uuidmap[item.uuid]
                if item.uuid in dist_matrix_uuids:
                    dist_matrix = dist_matrix_uuids[item.uuid]
                else:
                    dist_matrix = self.get_dist_matrix(system, max_cutoff)
                    dist_matrix_uuids[item.uuid] = dist_matrix

                equ_block, ref_block = self.make_energy_block(system, dist_matrix, self.modelcollection, self.external_repulsive_potentials)
                final_equ_block -= equ_block*item.coefficient
                final_ref_block -= ref_block*item.coefficient
                elec_reaction_energy -= system.elec_data['elec_energy']*item.coefficient


            final_ref_block += react_energy - elec_reaction_energy

            desc = EquationDescription('reaction', reaction_item, final_equ_block*reaction_item.weight, final_ref_block*reaction_item.weight)
            self.equation_descriptions.append(desc)

    def _build_additional_equations(self):
        for rep_name, items in sorted(self.evaluated_erepfit_input.equations['additional'].items(), key=lambda x:x[0]):
            shift = self.modelcollection.get_shift(rep_name)
            if (shift > -1):
                equ_block = np.zeros( (len(items), self.modelcollection.get_total_num_variables()) )
                ref_block = np.zeros( len(items) )
                for ind, item in enumerate(items):
                    r = utils.get_length_in_au(item)
                    row = self.modelcollection.get_row(rep_name, r, item.derivative)
                    lens = len(row)
                    equ_block[ind, shift:shift+lens] += row*item.weight
                    ref_block[ind] = item.value*item.weight

                desc = EquationDescription('additional', item, equ_block, ref_block)
                desc.rep_name = rep_name
                self.equation_descriptions.append(desc)
            else:
                continue

    def _build_force_equations(self, dist_matrix_uuids, max_cutoff):
        for item in self.evaluated_erepfit_input.equations['force']:
            if (item.weight < 1.0e-4):
                continue
            # Get the bond distances.
            system = self.uuidmap[item.uuid]
            if item.uuid in dist_matrix_uuids:
                dist_matrix = dist_matrix_uuids[item.uuid]
            else:
                dist_matrix = self.get_dist_matrix(system, max_cutoff)
                dist_matrix_uuids[item.uuid] = dist_matrix

            equ_block, indexes, ref_block = self.make_force_block(system, dist_matrix, self.modelcollection, self.external_repulsive_potentials)


            if (len(ref_block) == 0):
                continue

            total_forces = np.zeros( ref_block.shape )

            ref_e_forces = np.array([ system.elec_data['elec_forces'][i] for i in indexes ]).flatten()
            ref_block =  (total_forces -  ref_e_forces) - ref_block

            desc = EquationDescription('forces',item, equ_block*item.weight, ref_block*item.weight)
            self.equation_descriptions.append(desc)

    def _build_energy_equations(self, dist_matrix_uuids, max_cutoff):
        for item in self.evaluated_erepfit_input.equations['energy']:
            if (item.weight < 1.0e-4):
                continue
            # Get the bond distances.
            system = self.uuidmap[item.uuid]
            if item.uuid in dist_matrix_uuids:
                dist_matrix = dist_matrix_uuids[item.uuid]
            else:
                dist_matrix = self.get_dist_matrix(system, max_cutoff)
                dist_matrix_uuids[item.uuid] = dist_matrix

            e_atom = 0.0
            for line in system.geometry.coordinates:
                e_atom += self.evaluated_erepfit_input.atomic_energy[line[0]]

            e_ea_dft = utils.get_energy_in_au(item)
            e_rep_dftb = e_atom - system.elec_data['elec_energy'] - e_ea_dft

            equ_block, ref_block = self.make_energy_block(system, dist_matrix, self.modelcollection, self.external_repulsive_potentials)
            ref_block = -ref_block + e_rep_dftb

            desc = EquationDescription('energy', item, equ_block*item.weight, ref_block*item.weight)
            self.equation_descriptions.append(desc)

    def make_Q_R_S_matrix(self, dx, size):
        """ Get the continuity condition equations
        """
        mat_full = np.zeros( (size, size+1))
        mat_shift_full = np.zeros( (size, size+1))

        for i in range(0, size):
            for j in range(i, size+1):
                mat_shift_full[i, j] = 1.0*dx**(j-i)
                if (i==0):
                    mat_full[i, j] = 1.0
                else:
                    mat_full[i, j] = mat_full[i-1, j-1]*j

        final_mat = mat_full*mat_shift_full
        Q_mat = final_mat[0:size, 0:size]
        S_mat = final_mat[0:size, size-5:size+1]
        R_mat = np.diagflat(np.diag(Q_mat))*-1
        
        return Q_mat, R_mat, S_mat
    
    def make_Q_R_matrix(self, dx, nvar, order):
        nrows = order+1
        ncols = nvar
        mat_full = np.zeros( (nrows, ncols))
        mat_shift_full = np.zeros( (nrows, ncols))
        for i in range(0, nrows):
            for j in range(i, ncols):
                mat_shift_full[i, j] = 1.0*dx**(j-i)
                if (i==0):
                    mat_full[i, j] = 1.0
                else:
                    mat_full[i, j] = mat_full[i-1, j-1]*j

        final_mat = mat_full*mat_shift_full
        R_mat = np.diagflat(np.diag(final_mat))*-1
        
        return final_mat, R_mat

    # def solve_with_continuous_condition2(self, modelcollection, continuous_order=3):
        
    #     #avaliable only for Spline model
    #     if (not isinstance(modelcollection.rep_models[0], rep_spline.Spline4RepulsiveModel)):
    #         raise RuntimeError('The repulsive potential is not represented by Spline functions')
        
    #     #define the dimension of the matrices
    #     nsplines = sum([ len(x.spline_models) for x in modelcollection.rep_models])
    #     nvar = modelcollection.rep_models[0].spline_models[0].num_variables
        
    #     block_size = continuous_order+1
    #     w_block_size = nvar-block_size
    #     mat_C = np.zeros( ((block_size*nsplines), modelcollection.get_total_num_variables()) )
    #     vec_D = np.zeros( (block_size*nsplines) )               
        
    #     # print('Size of T and W ', mat_T.shape, mat_W.shape)
    #     # ref_vref = np.genfromtxt('../ref/vref.dat')
    #     # ref_eqmat = np.genfromtxt('../ref/eqmat.dat')

    #     ind = 0
    #     for spline_model in modelcollection.rep_models:
    #         knots = spline_model.knots
    #         for i in range(len(knots)-1):
    #             dx = knots[i+1] - knots[i]
    #             Q_mat, R_mat = self.make_Q_R_matrix(dx, nvar, continuous_order)
    #             mat_C[ind*block_size:(ind+1)*block_size, ind*nvar:(ind+1)*nvar] = Q_mat
                
    #             if (i < len(knots)-2):
    #                 mat_C[ind*block_size:(ind+1)*block_size, (ind+1)*nvar:(ind+1)*nvar+block_size] = R_mat

    #             ind += 1

    #     eqmat = np.concatenate([ x.equ_block for x in self.equation_descriptions ])
    #     vref = np.concatenate([ x.ref_block for x in self.equation_descriptions ])
        
    #     t, r, res, x, info = scipy.linalg.lapack.dgglse(eqmat, mat_C, vref, vec_D)

    #     print(x)
    #     modelcollection.resulting_repulsives = []
    #     shift = 0
    #     for ind in range(len(modelcollection.rep_names)):
    #         rep_name = modelcollection.rep_names[ind]
    #         rep_model = modelcollection.rep_models[ind]
    #         knots = rep_model.knots

    #         coeffs = []
    #         for i in range(len(knots)-1):
    #             coeffs.append( x[shift:shift+nvar] )
    #             shift += nvar
            
    #         spl_rep = spline_rep.RepulsivePotenial(rep_name, knots, coeffs)
    #         with open('{}.skf'.format(rep_name), 'w') as fout:
    #             print(spl_rep, file=fout)
    #         modelcollection.resulting_repulsives.append(spl_rep)

    def solve_with_continuous_condition(self, modelcollection, continuous_order=3):
        
        #avaliable only for Spline model
        if (not isinstance(modelcollection.rep_models[0], rep_spline.Spline4RepulsiveModel)):
            raise RuntimeError('The repulsive potential is not represented by Spline functions')

        #transform the equations using continuous conditions
        #define the dimension of the matrices
        nsplines = sum([ len(x.spline_models) for x in modelcollection.rep_models])
        nvar = modelcollection.rep_models[0].spline_models[0].num_variables
        block_size = continuous_order+1
        w_block_size = nvar-block_size
        mat_T = np.zeros( ((block_size*nsplines), (block_size*nsplines)))
        mat_W = np.zeros( ((block_size*nsplines), (w_block_size*nsplines)))

        ind = 0
        for spline_model in modelcollection.rep_models:
            knots = spline_model.knots
            for i in range(len(knots)-1):
                dx = knots[i+1] - knots[i]
                Q_mat, R_mat, S_mat = self.make_Q_R_S_matrix(dx, block_size)
                mat_T[ind*block_size:(ind+1)*block_size, ind*block_size:(ind+1)*block_size] = Q_mat
                if (i < len(knots)-2):
                    mat_T[ind*block_size:(ind+1)*block_size, (ind+1)*block_size:(ind+2)*block_size] = R_mat
                mat_W[ind*block_size:(ind+1)*block_size, (ind)*w_block_size:(ind+1)*w_block_size] = S_mat
                ind += 1

        eqmat = np.concatenate([ x.equ_block for x in self.equation_descriptions ])
        vref = np.concatenate([ x.ref_block for x in self.equation_descriptions ])
        
        #make eqmat separately for e3 and e4
        e3_index = sum([list(range(i, i+block_size )) for i in range(0, eqmat.shape[1]-block_size, nvar)], [])
        e3 = eqmat[:, e3_index]
        e4_index = sum([list(range(i, i+w_block_size )) for i in range(block_size, eqmat.shape[1], nvar)], [])
        e4 = eqmat[:, e4_index]       
        mat_T_inv_W = np.dot(np.linalg.inv(mat_T), -mat_W)
        mat_A  = np.dot(e3, mat_T_inv_W) + e4


        #solve with regularization solver
        reg = linear_model.LinearRegression()
        # reg = linear_model.Lasso(alpha = 0.01)
        reg.fit(mat_A, vref)

        #transform back for all coefficients
        a4 = reg.coef_
        # print(a4)
        a3 = np.dot(mat_T_inv_W, a4)


        modelcollection.resulting_repulsives = []
        vunknowns = []
        shift = 0
        for ind in range(len(modelcollection.rep_names)):
            rep_name = modelcollection.rep_names[ind]
            rep_model = modelcollection.rep_models[ind]
            knots = rep_model.knots

            coeffs = []
            for i in range(len(knots)-1):
                coeffs.append( np.append(a3[shift*block_size:shift*block_size+block_size], a4[shift:shift+w_block_size] ))
                vunknowns += a3[shift*block_size:shift*block_size+block_size].tolist()
                vunknowns += a4[shift:shift+w_block_size].tolist()
                shift += 1
            
            
            spl_rep = rep_spline.RepulsivePotenial(rep_name, knots, coeffs)
            modelcollection.resulting_repulsives.append(spl_rep)

        self.vunknowns = vunknowns
        

    def print_residuals(self):
        print("")
        print("Energy Equation Residuals")
        vunknowns = np.array(self.vunknowns)
        for desc in [x for x in self.equation_descriptions if x.equ_type == 'energy']:
            item = desc.equation_item
            res = (np.dot(desc.equ_block, vunknowns) - desc.ref_block) / item.weight
            print('{:10.4f} {:9s} {}'.format(utils.get_energy_from_au(res[0], item), item.unit, item.name))
        
        print("")
        print("Force Equation Residuals (Ha/bohr)")
        for desc in [x for x in self.equation_descriptions if x.equ_type == 'forces']:
            item = desc.equation_item
            res = (np.dot(desc.equ_block, vunknowns) - desc.ref_block) / item.weight
            print(item.name, item.weight)
            system = self.uuidmap[item.uuid]
            
            for i in range(0, int(len(res)/3)):
                print(' {:3s} {:20.12f} {:20.12f} {:20.12f}'.format(system.geometry.coordinates[i][0], 
                    res[i*3], res[i*3+1], res[i*3+2]))
            print("")

        print("")
        print("Addition Equation Residuals (Ha)")
        for desc in [x for x in self.equation_descriptions if x.equ_type == 'additional']:
            item = desc.equation_item
            rep_name = desc.rep_name
            res = np.dot(desc.equ_block, vunknowns) - desc.ref_block           
            for i in range(len(res)):
                print('  {:<6s} {:20.12f}'.format(rep_name, res[i]))

        print("")
        print("Reaction Equation Residuals (Ha)")
        for desc in [x for x in self.equation_descriptions if x.equ_type == 'reaction']:
            item = desc.equation_item
            res = (np.dot(desc.equ_block, vunknowns) - desc.ref_block ) / item.weight
            for i in range(len(res)):
                print('{:10.4f} {:9s} {}'.format(utils.get_energy_from_au(res[i], item), item.unit, item.name))
        print("")


    def solve(self, modelcollection, solver=linear_model.LinearRegression()):
        
        #avaliable only for Spline model
        if (not isinstance(modelcollection.rep_models[0], rep_basis.BasisRepulsiveModel)):
            raise RuntimeError('The repulsive potential is not represented by basis functions')
               
        eqmat = np.concatenate([ x.equ_block for x in self.equation_descriptions ])
        vref = np.concatenate([ x.ref_block for x in self.equation_descriptions ])

        solver.fit(eqmat, vref)
        # #transform back for all coefficients
        a4 = solver.coef_
        

        modelcollection.resulting_repulsives = []

        shift = 0
        for ind in range(len(modelcollection.rep_names)):
            rep_name = modelcollection.rep_names[ind]
            rep_model = modelcollection.rep_models[ind]
            nvar = rep_model.num_variables

            coeffs = a4[shift:shift+nvar] 
            shift += nvar

            newrep = rep_basis.PolynomialRepulsive(rep_name, rep_model.cutoff, coeffs, rep_model.minimal_order)
            spl_rep = newrep.to_spline4()
            modelcollection.resulting_repulsives.append(spl_rep)
        
        self.vunknowns = a4