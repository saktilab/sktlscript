#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy import stats
from collections import deque

def read_dump(filename):
    data = []
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                timestep = int(f.readline().strip())
                f.readline()  # ITEM: NUMBER OF ATOMS
                num_atoms = int(f.readline().strip())
                f.readline()  # ITEM: BOX BOUNDS
                f.readline()
                f.readline()
                f.readline()
                f.readline()  # ITEM: ATOMS id element xu yu zu
                atoms = []
                for _ in range(num_atoms):
                    atom_data = f.readline().split()
                    atom_id = int(atom_data[0])
                    element = atom_data[1]
                    x, y, z = map(float, atom_data[2:5])
                    atoms.append([atom_id, element, x, y, z])
                data.append((timestep, atoms))
    return data

def find_all_potential_hydroniums(atoms, max_oh_bond=1.5, min_oh_bond=0.8):
    o_atoms = np.array([(atom[0], *atom[2:5]) for atom in atoms if atom[1] == 'O'])
    h_atoms = np.array([atom[2:5] for atom in atoms if atom[1] == 'H'])

    if len(o_atoms) == 0 or len(h_atoms) == 0:
        return []

    distances = cdist(o_atoms[:, 1:], h_atoms)

    potential_hydroniums = []
    for i, o_data in enumerate(o_atoms):
        o_id, ox, oy, oz = o_data
        close_h = np.where((distances[i] <= max_oh_bond) & (distances[i] >= min_oh_bond))[0]
        if len(close_h) >= 2:  # Consider oxygen with at least 2 close hydrogens
            h_positions = h_atoms[close_h]
            potential_hydroniums.append({
                'O_id': int(o_id),
                'O_position': np.array([ox, oy, oz]),
                'H_positions': h_positions,
                'O_H_distances': distances[i][close_h],
                'num_close_h': len(close_h)
            })

    return potential_hydroniums

def calculate_refined_h_function(trajectory_data, max_oh_bond=1.5, min_oh_bond=0.8, persistence_check=5):
    h = [0]
    o_ids = []
    current_hydronium_id = None
    potential_transfer_queue = deque(maxlen=persistence_check)

    for timestep, atoms in trajectory_data:
        potential_hydroniums = find_all_potential_hydroniums(atoms, max_oh_bond, min_oh_bond)

        if not potential_hydroniums:
            if current_hydronium_id is not None:
                h.append(h[-1])
                o_ids.append(current_hydronium_id)
            continue

        # Sort potential hydroniums by number of close hydrogens, then by O-H distances
        potential_hydroniums.sort(key=lambda x: (-x['num_close_h'], np.min(x['O_H_distances'])))

        most_likely_hydronium = potential_hydroniums[0]

        if current_hydronium_id is None:
            current_hydronium_id = most_likely_hydronium['O_id']

        if most_likely_hydronium['O_id'] != current_hydronium_id:
            # Potential proton transfer
            potential_transfer_queue.append(most_likely_hydronium['O_id'])

            if len(potential_transfer_queue) == persistence_check and all(id == potential_transfer_queue[0] for id in potential_transfer_queue):
                # Confirmed proton transfer
                h.append(h[-1] + 1)
                current_hydronium_id = most_likely_hydronium['O_id']
                potential_transfer_queue.clear()
            else:
                h.append(h[-1])
        else:
            potential_transfer_queue.clear()
            h.append(h[-1])

        o_ids.append(current_hydronium_id)

    return np.array(h), np.array(o_ids)

def calculate_diffusion_coefficient(times, h, proton_transfer_length=2.5):
    # Ensure times and h have the same length
    min_length = min(len(times), len(h))
    times = times[:min_length]
    h = h[:min_length]

    # Calculate the gradient of h
    gradient_h = np.gradient(h, times)

    # Calculate mean squared gradient
    mean_squared_gradient = np.mean(gradient_h**2)

    # Calculate diffusion coefficient
    D = (mean_squared_gradient * proton_transfer_length**2) / 6

    # Perform linear regression for visualization purposes
    slope, intercept, r_value, _, std_err = stats.linregress(times, h)

    return D, slope, intercept, r_value**2, std_err

def main(args):
    # Read trajectory
    print("Reading trajectory...")
    data = read_dump(args.dump_file)

    # Calculate refined h function and get oxygen IDs
    print("Calculating refined h function...")
    h, o_ids = calculate_refined_h_function(data, max_oh_bond=args.max_oh_bond, min_oh_bond=args.min_oh_bond, persistence_check=args.persistence_check)
    times = np.arange(len(h)) * args.timestep / 1000  # Convert to ps

    # Ensure times, h, and o_ids have the same length
    min_length = min(len(times), len(h), len(o_ids))
    times = times[:min_length]
    h = h[:min_length]
    o_ids = o_ids[:min_length]

    # Calculate diffusion coefficient
    D, slope, intercept, r_squared, std_err = calculate_diffusion_coefficient(times, h, args.proton_transfer_length)

    print(f"Diffusion coefficient: {D:.6f} Å²/ps")
    print(f"R-squared value of linear fit: {r_squared:.6f}")

    # Plot h function with linear fit
    print("Generating plots...")
    plt.figure(figsize=(12, 12))

    plt.subplot(311)
    plt.plot(times, h, 'b.', label='h(t)')
    plt.plot(times, slope * times + intercept, 'r-', label='Linear fit')
    plt.xlabel('Time (ps)')
    plt.ylabel('h(t)')
    plt.title('Hydronium Ion Parameter h(t) with Linear Fit')
    plt.legend()
    plt.grid(True)

    # Plot gradient of h function
    plt.subplot(312)
    gradient = np.gradient(h, times)
    plt.plot(times, gradient, 'g-')
    plt.xlabel('Time (ps)')
    plt.ylabel('dh/dt')
    plt.title('Gradient of h(t)')
    plt.grid(True)

    # Plot oxygen ID
    plt.subplot(313)
    plt.plot(times, o_ids, 'r-')
    plt.xlabel('Time (ps)')
    plt.ylabel('Oxygen Atom ID')
    plt.title('Oxygen Atom ID of Hydronium Ion')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(args.output_plot)
    plt.close()

    print(f"Plots saved as {args.output_plot}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze hydronium ion shuttling from LAMMPS trajectory")
    parser.add_argument("dump_file", type=str, help="Path to the LAMMPS dump file")
    parser.add_argument("--timestep", type=float, default=1.0, help="Simulation timestep in fs")
    parser.add_argument("--max_oh_bond", type=float, default=1.5, help="Maximum O-H bond length to consider (Angstrom)")
    parser.add_argument("--min_oh_bond", type=float, default=0.8, help="Minimum O-H bond length to consider (Angstrom)")
    parser.add_argument("--proton_transfer_length", type=float, default=2.5, help="Average proton transfer length (Angstrom)")
    parser.add_argument("--persistence_check", type=int, default=5, help="Number of consecutive timesteps to confirm a proton transfer")
    parser.add_argument("--output_plot", type=str, default="hydronium_analysis_plot.png", help="Output file name for the plots")

    args = parser.parse_args()
    main(args)
