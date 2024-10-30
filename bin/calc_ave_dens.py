#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def read_lammps_log(log_file):
    densities = []
    steps = []
    with open(log_file, 'r') as f:
        for line in f:
            if line.startswith('Step'):
                headers = line.split()
                density_index = headers.index('Density')
                step_index = headers.index('Step')
            elif 'Loop' in line:
                break
            elif line[0].isdigit():
                values = line.split()
                densities.append(float(values[density_index]))
                steps.append(int(values[step_index]))
    return np.array(steps), np.array(densities)

def detect_convergence(steps, densities,threshold, window=10000):
    for i in range(len(densities) - window):
        slope, _, _, _, _ = linregress(steps[i:i+window], densities[i:i+window])
        if abs(slope) < threshold:
            return i, np.mean(densities[i:])
    return -1, None

def find_closest_frame(steps, densities, target_density):
    closest_index = np.argmin(np.abs(densities - target_density))
    return steps[closest_index], densities[closest_index]

def extract_frame_from_dump(dump_file, target_step):
    frame_data = []
    reading_frame = False
    with open(dump_file, 'r') as f:
        for line in f:
            if 'ITEM: TIMESTEP' in line:
                step = int(next(f).strip())
                if step == target_step:
                    reading_frame = True
                elif reading_frame:
                    break
            if reading_frame:
                frame_data.append(line)
    return frame_data

def convert_dump_to_data(frame_data):
    data_file = []
    atom_data = []
    box_bounds = []
    num_atoms = 0

    for line in frame_data:
        if 'ITEM: NUMBER OF ATOMS' in line:
            num_atoms = int(next(frame_data))
        elif 'ITEM: BOX BOUNDS' in line:
            box_bounds = [next(frame_data).split() for _ in range(3)]
        elif 'ITEM: ATOMS' in line:
            headers = line.split()[2:]
            for _ in range(num_atoms):
                atom_data.append(next(frame_data).split())

    # Write LAMMPS data file header
    data_file.append("LAMMPS data file from dump frame\n")
    data_file.append(f"{num_atoms} atoms\n")
    data_file.append(f"{len(set(atom[1] for atom in atom_data))} atom types\n")

    # Write box dimensions
    for i, (lo, hi) in enumerate(box_bounds):
        data_file.append(f"{lo} {hi} xlo xhi\n" if i == 0 else
                         f"{lo} {hi} ylo yhi\n" if i == 1 else
                         f"{lo} {hi} zlo zhi\n")

    # Write atom data
    data_file.append("\nAtoms\n\n")
    for atom in atom_data:
        data_file.append(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]} {atom[4]}\n")

    return data_file

def main(log_file, dump_file,threshold):
    steps, densities = read_lammps_log(log_file)

    convergence_index, avg_density = detect_convergence(steps, densities,threshold)

    if convergence_index != -1:
        print(f"Density converged at step {steps[convergence_index]}")
        print(f"Average density after convergence: {avg_density:.6f}")

        closest_step, closest_density = find_closest_frame(steps, densities, avg_density)
        print(f"Frame with closest density: Step {closest_step}, Density {closest_density:.6f}")

        frame_data = extract_frame_from_dump(dump_file, closest_step)
        if frame_data:
            data_file_content = convert_dump_to_data(frame_data)
            with open('closest_frame.data', 'w') as f:
                f.writelines(data_file_content)
            print("Frame with closest density converted and saved as 'closest_frame.data'")
        else:
            avg_density = np.mean(densities[-1000:])  # Average over last 1000 steps
            print(f"Average density over last 1000 steps: {avg_density:.6f}")
    # ... proceed with finding closest frame and outputting data file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze LAMMPS log and dump files')
    parser.add_argument('-l','--log_file', type=str, help='Path to the LAMMPS log file')
    parser.add_argument('-d','--dump_file', type=str, help='Path to the LAMMPS dump file')
    parser.add_argument('-t','--threshold', type=float, help='Convergence threshold')
    args = parser.parse_args()

    main(args.log_file, args.dump_file,args.threshold)
