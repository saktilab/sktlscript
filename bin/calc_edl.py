#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define system parameters
eps0 = 8.854e-12  # permittivity of free space in F/m
e_charge = 1.602e-19  # elementary charge in C
k_B = 1.381e-23  # Boltzmann constant in J/K
T = 300  # temperature in K
ion_conc = 0.1  # ionic concentration in M
sigma = 0.1  # surface charge density in C/m^2

# Load XYZ trajectory data
# This assumes that the XYZ file has three columns for the x, y, and z coordinates, respectively
positions = np.loadtxt('trajectory.xyz', skiprows=2, usecols=(1, 2, 3))

# Define grid points for potential calculation
x = np.linspace(-10, 10, 101)  # grid points in x direction
y = np.linspace(-10, 10, 101)  # grid points in y direction
z = np.linspace(-10, 10, 101)  # grid points in z direction
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')  # create 3D grid

# Calculate Coulomb potential at each grid point
r = np.sqrt((xx - positions[:, 0][:, np.newaxis, np.newaxis])**2 + (yy - positions[:, 1][:, np.newaxis, np.newaxis])**2 + (zz - positions[:, 2][:, np.newaxis, np.newaxis])**2)  # distance from each grid point to each atom
q = e_charge * ion_conc  # charge density in solution
potential = 1 / (4 * np.pi * eps0 * r) * (np.exp(-e_charge * sigma / (k_B * T)) - np.exp(e_charge * sigma / (k_B * T))) / (np.exp(-e_charge * q / (eps0 * k_B * T)) + np.exp(e_charge * q / (eps0 * k_B * T)))
total_potential = np.sum(potential, axis=0)  # sum up the potential from all atoms

# Plot the potential surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, facecolors=plt.cm.viridis(total_potential), shade=False)
plt.show()
