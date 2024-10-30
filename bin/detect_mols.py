#!/usr/bin/env python
import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use the 'Agg' backend
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import sys
from matplotlib.backends.backend_pdf import PdfPages

# Average van der Waals radius (in Angstroms)
average_vdw_radius = 1.5  # This is an approximate value, adjust as needed

def compact_formula(formula):
    elements = {}
    for element in formula:
        elements[element] = elements.get(element, 0) + 1
    return ''.join(f"{e}{n if n > 1 else ''}" for e, n in elements.items())

def detect_molecules(positions, elements, distance_cutoff):
    distances = squareform(pdist(positions))
    adjacency_matrix = distances < distance_cutoff
    n = len(positions)
    visited = np.zeros(n, dtype=bool)
    molecules = []

    for i in range(n):
        if not visited[i]:
            molecule = []
            stack = [i]
            while stack:
                atom = stack.pop()
                if not visited[atom]:
                    visited[atom] = True
                    molecule.append(atom)
                    stack.extend(np.where(adjacency_matrix[atom] & ~visited)[0])
            molecules.append(molecule)

    chemical_formulas = [compact_formula(''.join(sorted(elements[molecule]))) for molecule in molecules]
    return chemical_formulas

# Load the XYZ trajectory
u = mda.Universe(sys.argv[1])

# Dictionary to store the time course of each species
species_time_course = defaultdict(list)

# Get all elements
elements = np.array([atom.name for atom in u.atoms])

# Time step in picoseconds (10 fs = 0.01 ps)
time_step = 0.01

# Set distance cutoff (2 * average_vdw_radius)
distance_cutoff = float(sys.argv[2]) 

# Iterate through each frame with a progress bar
for frame, ts in enumerate(tqdm(u.trajectory, desc="Processing frames", unit="frame")):
    positions = u.atoms.positions
    molecules = detect_molecules(positions, elements, distance_cutoff)
    unique, counts = np.unique(molecules, return_counts=True)
    molecule_counts = dict(zip(unique, counts))


    for species in set(species_time_course.keys()) | set(molecule_counts.keys()):
        species_time_course[species].append(molecule_counts.get(species, 0))

# Convert time course data to numpy arrays
for species in species_time_course:
    species_time_course[species] = np.array(species_time_course[species])

# Create time array in picoseconds
time = np.arange(len(u.trajectory)) * time_step

# Ensure all arrays have the same length as the time array
for species in species_time_course:
    if len(species_time_course[species]) < len(time):
        species_time_course[species] = np.pad(species_time_course[species],
                                              (0, len(time) - len(species_time_course[species])),
                                              'constant', constant_values=0)
    elif len(species_time_course[species]) > len(time):
        species_time_course[species] = species_time_course[species][:len(time)]

# Calculate the variance of each species
species_variance = {species: np.var(counts) for species, counts in species_time_course.items()}

# Select the most important species 
important_species = sorted(species_variance, key=species_variance.get, reverse=True)[:4]
important_species = list(dict.fromkeys(important_species))[:6]  # Remove duplicates and limit to 6

# Create the plot and save as PDF
with PdfPages('species_time_course.pdf') as pdf:
    plt.figure(figsize=(12, 7))  # Increased figure size
    for species in important_species:
        if species in species_time_course:
            plt.plot(time, species_time_course[species], label=species, linewidth=1)

    plt.xlabel('Time [ps]', fontsize=14)
    plt.ylabel('Number of molecules', fontsize=14)
    leg = plt.legend(fontsize=14, loc='center left', bbox_to_anchor=(1, 0.5))
    
    for line in leg.get_lines():
        line.set_linewidth(1.5)

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=8, direction='in')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

# Write summary statistics to a file
with open('molecular_species_statistics.txt', 'w') as f:
    f.write("Summary Statistics:\n")
    for species, counts in species_time_course.items():
        f.write(f"\n{species}:\n")
        f.write(f"  Average count: {np.mean(counts):.2f}\n")
        f.write(f"  Maximum count: {np.max(counts)}\n")
        f.write(f"  Minimum count: {np.min(counts)}\n")
        f.write(f"  Final count: {counts[-1]}\n")
        f.write(f"  Variance: {np.var(counts):.2f}\n")
        f.write(f"  Data points: {len(counts)}\n")

print("Done!")
