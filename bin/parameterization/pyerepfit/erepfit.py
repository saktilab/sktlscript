#!/usr/bin/env python3

import argparse
import os
import sys
import json
import jsonschema
import yaml
import jsonpickle

from inputdata import ErepfitInput
import evaluator
from rep_basis import PolynomialRepulsive, BasisRepulsiveModel
from rep_spline import RepulsivePotenial, Spline4RepulsiveModel
import rep_model
import logging

VERSION = 0.1

import os

def which(file):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
            return os.path.join(path, file)

    return file

def print_header():
    print(
"""
===============================================================================

    PyErepfit v{}
    Powered by Chien-Pin Chou

===============================================================================
""".format(VERSION)
    )

            
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", dest='yaml', help='Parse input as YAML format', action='store_true')
    parser.add_argument("-e", dest='evaluate_only', help='Evaluate the DFTB only', action='store_true')
    parser.add_argument("-f", dest='fitting_only', help='Fitting the repulsive only', action='store_true')
    parser.add_argument("input_file_path", help="Input file path")


    opt = parser.parse_args(sys.argv[1:])

    print_header()

    print('Reading input schema.')
    with open(os.path.join(os.path.dirname(__file__), 'schema', 'erepfit_input_schema.json'), 'r') as f:
        input_schema = json.load(f)

    if (opt.yaml):
        print('Reading input file as YAML')
        with open(opt.input_file_path, 'r') as f:
            input_file = yaml.load(f)
    else:
        print('Reading input file as JSON')
        with open(opt.input_file_path, 'r') as f:
            input_file = json.load(f)

    print('Validate input file.')
    jsonschema.validate(input_file, input_schema)

    erepfit_input = ErepfitInput(input_file)

    erepfit_input.electronic_slater_koster_files.adjust_path(os.path.abspath(os.curdir))


    rep_models = []
    rep_names = []
    poten_type = erepfit_input.options.get("potential_type")
    
    for k, v in erepfit_input.potential_grids.items():
        if (k in erepfit_input.external_repulsive_potentials):
            logging.info('Repulsive {} found in external repulsuive potentials, will not be fitted.'.format(k))
            continue
        if (poten_type == "polynomial"):
            cutoff = v['knots'][-1]
            model = BasisRepulsiveModel(cutoff, minimal_order=3, order=9)
        else: 
            model = Spline4RepulsiveModel(v['knots'])
        rep_names.append(k)
        rep_models.append(model)

    collection = rep_model.RepulsiveModelCollection(rep_names, rep_models)
    print("")
    print(collection.description())

    if (opt.evaluate_only or not opt.fitting_only):
        print('Testing DFTB implementation for parameterization...')
        dftb_path = erepfit_input.options['toolchain']['path']
        if not os.path.isabs(dftb_path):
            dftb_path = os.path.expanduser(dftb_path)
        dftb_path = os.path.expandvars(dftb_path)
        if not os.path.isabs(dftb_path):
            dftb_path = which(dftb_path)

        succ = evaluator.test_dftb_binary(dftb_path)
        if (succ):
            print('  {} => working'.format(dftb_path))
        else:
            print('  {} => not working'.format(dftb_path))
            sys.exit(1)
        
        
        evaluator.evaluate(erepfit_input)
        print('\nEvaluate DFTB...done')

        jsonpickle.set_encoder_options('simplejson', indent=4)
        jsonpickle.set_encoder_options('json', indent=4)
        
        with open('evaluated.json', 'w') as f:
            js = jsonpickle.encode(erepfit_input, unpicklable=False)
            print(js, file=f)
    

    if (opt.fitting_only or not opt.evaluate_only):
        solver = rep_model.RepulsivePotentialEquationBuilder(erepfit_input)
        solver.build_linear_equations(collection)
        
        if (poten_type == "polynomial"):
            solver.solve(collection)
        else:
            solver.solve_with_continuous_condition(collection, 3)
        
        solver.print_residuals()

        print('Resulting Repulsive Potentials:')
        output_prefix = erepfit_input.options.get("output_prefix", "./")
        os.makedirs(output_prefix, exist_ok=True)
        for ind in range(len(collection.rep_names)):
            rep_name = collection.rep_names[ind]
            spl_rep = collection.resulting_repulsives[ind]
            reversed_rep_name = '{}-{}'.format(*reversed(rep_name.split('-')))

            with open(os.path.join(output_prefix, '{}.skf'.format(rep_name)), 'w') as fout:
                print(spl_rep, file=fout)
            if (reversed_rep_name != rep_name):
                with open(os.path.join(output_prefix, '{}.skf'.format(reversed_rep_name)), 'w') as fout:
                    print(spl_rep, file=fout)
            
            print(rep_name)
            print(spl_rep)

        print('End of pyerepfit Program')        
    
            




