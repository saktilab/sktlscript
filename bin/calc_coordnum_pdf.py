#!/usr/bin/env python
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import sys
def read_xyz(filename):
    """Read XYZ file and return atom positions and types."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    frames = []
    i = 0
    while i < len(lines):
        n_atoms = int(lines[i])
        i += 2  # Skip the comment line
        frame = []
        for _ in range(n_atoms):
            atom_data = lines[i].split()
            atom_type = atom_data[0]
            position = [float(x) for x in atom_data[1:4]]
            frame.append((atom_type, position))
            i += 1
        frames.append(frame)
    return frames

def calculate_coordination_numbers(frame, central_elements, neighbor_elements, cutoff):
    """Calculate coordination numbers for specified elements."""
    positions = np.array([atom[1] for atom in frame])
    types = [atom[0] for atom in frame]
    tree = cKDTree(positions)

    coordination_numbers = []
    for i, (atom_type, pos) in enumerate(frame):
        if atom_type in central_elements:
            neighbors = tree.query_ball_point(pos, cutoff)
            coord_num = sum(1 for j in neighbors if i != j and types[j] in neighbor_elements)
            coordination_numbers.append(coord_num)

    return coordination_numbers

def coordination_histogram(trajectory, central_elements, neighbor_elements, cutoff, max_coord=12):
    """Calculate histogram of coordination numbers for specified elements."""
    all_coordination_numbers = []
    for frame in trajectory:
        coordination_numbers = calculate_coordination_numbers(frame, central_elements, neighbor_elements, cutoff)
        all_coordination_numbers.extend(coordination_numbers)

    hist, bin_edges = np.histogram(all_coordination_numbers, bins=range(max_coord + 2), density=True)
    return hist*100, bin_edges

def plot_histogram(hist, bin_edges, central_elements, neighbor_elements):
    """Plot the histogram of coordination numbers."""
    plt.figure(figsize=(10, 6))
    plt.bar(bin_edges[:-1], hist, width=1, edgecolor='black',color='blue')
    plt.xlabel('Coordination Number')
    plt.ylabel('Probability [%]')
    #plt.title(f'Distribution of Coordination Numbers\n{", ".join(central_elements)} coordinated by {", ".join(neighbor_elements)}')
    plt.xticks(range(len(hist)))
    #plt.grid(axis='y', alpha=0.3)
    plt.savefig('coordination_histogram.png',dpi=1000)
    plt.close()

if __name__ == "__main__":
    filename = sys.argv[1]
    cutoff = float(sys.argv[4])  # Adjust this cutoff distance as needed for your system
    central_elements = [sys.argv[2]]  # Elements to analyze as central atoms
    neighbor_elements = [sys.argv[3]]  # Elements to consider as neighbors

    # Read trajectory
    trajectory = read_xyz(filename)
    print(f"Read {len(trajectory)} frames from the trajectory.")

    # Calculate and plot histogram
    hist, bin_edges = coordination_histogram(trajectory, central_elements, neighbor_elements, cutoff)
    plot_histogram(hist, bin_edges, central_elements, neighbor_elements)

    print("Coordination number histogram calculation complete.")
    print("Plot saved as 'coordination_histogram.png'.")

    # Print the histogram data
    for i, prob in enumerate(hist):
        print(f"Coordination number {i}: {prob:.4f}")
