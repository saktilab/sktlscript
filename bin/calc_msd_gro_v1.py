#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy import stats

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

def find_hydronium(atoms, max_oh_bond=1.3, min_oh_bond=0.8):
    o_atoms = np.array([(atom[0], *atom[2:5]) for atom in atoms if atom[1] == 'O'])
    h_atoms = np.array([atom[2:5] for atom in atoms if atom[1] == 'H'])

    if len(o_atoms) == 0 or len(h_atoms) == 0:
        return None

    distances = cdist(o_atoms[:, 1:], h_atoms)

    for i, o_data in enumerate(o_atoms):
        o_id, ox, oy, oz = o_data
        close_h = np.where((distances[i] <= max_oh_bond) & (distances[i] >= min_oh_bond))[0]
        if len(close_h) == 3:
            h_positions = h_atoms[close_h]
            return {
                'O_id': int(o_id),
                'O_position': np.array([ox, oy, oz]),
                'H_positions': h_positions,
                'O_H_distances': distances[i][close_h],
                'center_of_mass': np.mean(np.vstack(([ox, oy, oz], h_positions)), axis=0)
            }

    return None

def calculate_h_function(hydronium_data):
    h = [0]
    o_ids = [hydronium_data[0]['O_id']]

    for i in range(1, len(hydronium_data)):
        if hydronium_data[i]['O_id'] != hydronium_data[i-1]['O_id']:
            # Shuttling event occurred
            delta_h = 1 if hydronium_data[i]['O_id'] > hydronium_data[i-1]['O_id'] else -1
            h.append(h[-1] + delta_h)
        else:
            h.append(h[-1])
        o_ids.append(hydronium_data[i]['O_id'])

    return np.array(h), np.array(o_ids)

def calculate_diffusion_coefficient(times, h):
    # Perform linear regression
    slope, _, r_value, _, std_err = stats.linregress(times, h)

    # The diffusion coefficient is the slope divided by 6 times square of O-O distance
    D = slope/6*2.5**2

    return D, slope, r_value**2, std_err

def main(args):
    # Read trajectory
    print("Reading trajectory...")
    data = read_dump(args.dump_file)

    # Find hydronium positions
    print("Identifying hydronium ions...")
    hydronium_data = []
    prev_hydronium = None
    for timestep, atoms in data:
        hydronium = find_hydronium(atoms, max_oh_bond=args.max_oh_bond, min_oh_bond=args.min_oh_bond)
        if hydronium is not None:
            hydronium_data.append(hydronium)
            prev_hydronium = hydronium
        else:
            if prev_hydronium is None:
                print(f"Warning: No hydronium found at timestep {timestep} and no previous hydronium data available.")
                continue
            print(f"Warning: No hydronium found at timestep {timestep}. Using previous hydronium data.")
            hydronium_data.append(prev_hydronium)

    if len(hydronium_data) == 0:
        print("Error: No hydronium ions found in the trajectory.")
        return

    # Calculate h function and get oxygen IDs
    print("Calculating h function...")
    h, o_ids = calculate_h_function(hydronium_data)
    times = np.arange(len(h)) * args.timestep / 1000  # Convert to ps

    # Calculate diffusion coefficient using linear regression
    D, slope, r_squared, std_err = calculate_diffusion_coefficient(times, h)

    print(f"Diffusion coefficient: {D:.6f} ± {std_err:.6f} Å²/ps")
    print(f"R-squared value: {r_squared:.6f}")

    # Plot h function with regression line
    print("Generating plots...")
    plt.figure(figsize=(12, 12))

    plt.subplot(311)
    plt.plot(times, h, 'b.', label='h(t)')
    plt.plot(times, slope * times, 'r-', label='Regression line')
    plt.xlabel('Time (ps)')
    plt.ylabel('h(t)')
    plt.title('Hydronium Ion Parameter h(t) with Regression Line')
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
    parser.add_argument("--max_oh_bond", type=float, default=1.3, help="Maximum O-H bond length to consider (Angstrom)")
    parser.add_argument("--min_oh_bond", type=float, default=0.8, help="Minimum O-H bond length to consider (Angstrom)")
    parser.add_argument("--output_plot", type=str, default="hydronium_analysis_plot.png", help="Output file name for the plots")

    args = parser.parse_args()
    main(args)
