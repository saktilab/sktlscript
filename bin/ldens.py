#!/usr/bin/env python3
import sys
import matplotlib
# matplotlib.use('PS')
# matplotlib.rc('text', usetex=True)
matplotlib.use('Qt5agg')
import matplotlib.pyplot as plt
import MDAnalysis
import MDAnalysis.analysis.rdf
import numpy as np
import pymatgen as mg
import argparse
import importlib
import builtins
import math
from lmfit import Model


parser = argparse.ArgumentParser(description='Linear Density')
parser.add_argument('-s', '--start', default=1, type=int)
parser.add_argument('-e', '--endstep', default=-1, type=int)
parser.add_argument('-x', '--axis', default='z', choices=['x', 'y', 'z'])

parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')

parser.add_argument('filename', metavar='filename', type=str)
parser.add_argument('-o', '--output', type=str)

parser.add_argument('-d', '--resolution', default=0.25, type=float, help='Resolution for bins')

opt = parser.parse_args(sys.argv[1:])



# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Linear Density --- :mod:`MDAnalysis.analysis.lineardensity`
===========================================================

A tool to compute mass and charge density profiles along the three
cartesian axes of the simulation cell. Works only for orthorombic,
fixed volume cells (thus for simulations in canonical NVT ensemble).
"""

import os.path as path

import numpy as np

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.util import deprecate

class LinearDensity(AnalysisBase):
    """Linear density profile

    Parameters
    ----------
    selection : AtomGroup
          any atomgroup
    grouping : str {'atoms', 'residues', 'segments', 'fragments'}
          Density profiles will be computed on the center of geometry
          of a selected group of atoms ['atoms']
    binsize : float
          Bin width in Angstrom used to build linear density
          histograms. Defines the resolution of the resulting density
          profile (smaller --> higher resolution) [0.25]
    start : int
          The frame to start at [0]
    stop : int
          The frame to end at [-1]
    step : int
          The step size through the trajectory in frames [0]
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``; the
          default is ``False``.

    Attributes
    ----------
    results : dict
          Keys 'x', 'y', and 'z' for the three directions. Under these
          keys, find 'pos', 'pos_std' (mass-weighted density and
          standard deviation), 'char', 'char_std' (charge density and
          its standard deviation), 'slice_volume' (volume of bin).

    Example
    -------
    First create a LinearDensity object by supplying a selection,
    then use the :meth:`run` method::

      ldens = LinearDensity(selection)
      ldens.run()


    .. versionadded:: 0.14.0

    """

    def __init__(self, selection, grouping='atoms', binsize=0.25, **kwargs):
        super(LinearDensity, self).__init__(selection.universe.trajectory,
                                            **kwargs)
        # allows use of run(parallel=True)
        self._ags = [selection]
        self._universe = selection.universe

        self.binsize = binsize

        # group of atoms on which to compute the COM (same as used in
        # AtomGroup.wrap())
        self.grouping = grouping

        # Dictionary containing results
        self.results = {'x': {'dim': 0}, 'y': {'dim': 1}, 'z': {'dim': 2}}
        # Box sides
        self.dimensions = self._universe.dimensions[:3]
        self.volume = np.prod(self.dimensions)
        # number of bins
        bins = (self.dimensions // self.binsize).astype(int) + 1

        # Here we choose a number of bins of the largest cell side so that
        # x, y and z values can use the same "coord" column in the output file
        self.nbins = bins.max()
        slices_vol = self.volume / bins

        self.keys = ['pos', 'pos_std']

        # Initialize results array with zeros
        for dim in self.results:
            idx = self.results[dim]['dim']
            self.results[dim].update({'slice volume': slices_vol[idx]})
            for key in self.keys:
                self.results[dim].update({key: np.zeros(self.nbins)})

        # Variables later defined in _prepare() method
        self.masses = None
        self.totalmass = None

    def _prepare(self):
        # group must be a local variable, otherwise there will be
        # issues with parallelization
        group = getattr(self._ags[0], self.grouping)

        # Get masses and charges for the selection
        try:  # in case it's not an atom
            self.masses = np.array([elem.total_mass() for elem in group])
        except AttributeError:  # much much faster for atoms
            self.masses = self._ags[0].masses

        self.totalmass = np.sum(self.masses)

    def _single_frame(self):
        self.group = getattr(self._ags[0], self.grouping)
        self._ags[0].wrap(compound=self.grouping)

        # Find position of atom/group of atoms
        if self.grouping == 'atoms':
            positions = self._ags[0].positions  # faster for atoms
        else:
            # COM for res/frag/etc
            positions = np.array([elem.centroid() for elem in self.group])

        for dim in ['x', 'y', 'z']:
            idx = self.results[dim]['dim']

            key = 'pos'
            key_std = 'pos_std'
            # histogram for positions weighted on masses
            hist, _ = np.histogram(positions[:, idx],
                                   weights=self.masses,
                                   bins=self.nbins,
                                   range=(0.0, max(self.dimensions)))

            self.results[dim][key] += hist
            self.results[dim][key_std] += np.square(hist)


    def _conclude(self):
        k = 6.02214076e-1  # divide by avodagro and convert from A3 to cm3

        # Average results over the  number of configurations
        for dim in ['x', 'y', 'z']:
            for key in ['pos', 'pos_std']:
                self.results[dim][key] /= self.n_frames
            # Compute standard deviation for the error
            self.results[dim]['pos_std'] = np.sqrt(self.results[dim][
                'pos_std'] - np.square(self.results[dim]['pos']))
            

        for dim in ['x', 'y', 'z']:
            self.results[dim]['pos'] /= self.results[dim]['slice volume'] * k
            self.results[dim]['pos_std'] /= self.results[dim]['slice volume'] * k

    def _add_other_results(self, other):
        # For parallel analysis
        results = self.results
        for dim in ['x', 'y', 'z']:
            key = 'pos'
            key_std = 'pos_std'
            results[dim][key] += other[dim][key]
            results[dim][key_std] += other[dim][key_std]


outfile = sys.stdout
make_plot = False
if (opt.output):
    outfile = open('{}.out'.format(opt.output), 'w')
    make_plot = True

lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)
(lx, ly, lz), (a, b, c) = box.lengths_and_angles


u = MDAnalysis.Universe(opt.filename, format='XYZ')
u.dimensions = np.array([lx, ly, lz, a, b, c], dtype=np.float32)


ag1 = u.select_atoms('all')

ldens = LinearDensity(ag1, verbose=True, binsize=opt.resolution)
ldens.run(start=opt.start, stop=opt.endstep)


box_len = ldens.dimensions["xyz".index(opt.axis)]
bins = np.linspace(0.0, box_len, num=ldens.nbins)

output = [bins]

output.append(ldens.results[opt.axis]['pos'])
output.append(ldens.results[opt.axis]['pos_std'])

density = ldens.totalmass / ldens.volume

final_nbins = len(bins) // 2 
final_bins = bins[final_nbins:len(bins)]-bins[final_nbins]

final_values = ldens.results[opt.axis]['pos'][final_nbins:len(bins)]

for i in range(final_nbins):
    ind = final_nbins - i 
    final_values[i] += ldens.results[opt.axis]['pos'][ind]

for i in range(final_nbins):
    final_values[i] /= 2.0

def func(x, rho, z_gds, delta, c1, c2):
    return [(rho/2.0)*(1.0-math.tanh(c1+c2*(z-z_gds)/delta)) for z in x]

gmodel = Model(func)

params = gmodel.make_params(rho=1.0, z_gds=10.0, delta=1.0, c1=0.0, c2=1.0)
result = gmodel.fit(final_values, params, x=final_bins.tolist())
print(f'# Fitted Density: {result.values["rho"]}')
print(result.fit_report(), file=sys.stderr)

#np.savetxt(outfile, np.column_stack(output), fmt='%10.5f')

final_xs =final_bins.tolist()
rho, z_gds, delta, c1, c2 = result.values["rho"], result.values["z_gds"], result.values["delta"], result.values["c1"], result.values["c2"]
y = func(final_xs, rho, z_gds, delta, c1, c2) 
for r, calc_y, ref_y in zip(final_xs, y, final_values):
    print('{:12.5f} {:12.5f} {:12.5f}'.format(r, calc_y, ref_y))
