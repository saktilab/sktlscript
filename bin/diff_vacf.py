#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Define constants and parameters
kB = 1.38e-23  # Boltzmann constant
T = 300  # temperature in Kelvin
eta = 0.001  # viscosity of the solvent in Pa*s
ion_radius = 1e-9  # radius of the ion in meters
time_step = 1e-15  # time step in seconds

# Load position file
with open('trajectory.xyz', 'r') as f:
    num_ions = int(f.readline().strip())
    positions = np.zeros((num_ions, 3))
    for i in range(num_ions):
        line = f.readline().strip().split()
        positions[i] = [float(line[1]), float(line[2]), float(line[3])]

# Load velocity file
with open('velocities.xyz', 'r') as f:
    velocities = np.zeros((num_ions, 3))
    for i in range(num_ions):
        line = f.readline().strip().split()
        velocities[i] = [float(line[1]), float(line[2]), float(line[3])]

# Calculate center-of-mass velocity for each ion
com_velocities = np.mean(velocities, axis=0)

# Calculate velocity autocorrelation function
num_steps = len(positions)
vacf = np.zeros(num_steps)
for t in range(num_steps):
    for tau in range(num_steps - t):
        vacf[tau] += np.dot(velocities[t], velocities[t+tau])
    vacf[tau] /= (num_steps - t)

# Normalize velocity autocorrelation function
vacf /= np.amax(vacf)

# Calculate diffusion coefficient
D = kB*T/(6*np.pi*eta*ion_radius) * np.sum(vacf)*time_step

# Plot velocity autocorrelation function
time_range = np.arange(0, num_steps*time_step, time_step)
plt.plot(time_range[:len(vacf)], vacf)
plt.xlabel('Time (s)')
plt.ylabel('Velocity autocorrelation function')
plt.show()

print('Diffusion coefficient: {:.2e} m^2/s'.format(D))
