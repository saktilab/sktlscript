#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import io

def plot_lammps_log(log_file,dt):
    # Read the log file, skipping lines until the header is found
    with open(log_file, 'r') as f:
        lines = f.readlines()
        start_line = next(i for i, line in enumerate(lines) if 'CPU' in line and 'Step' in line)
        # Collect data until 'Loop' is encountered
        data_lines = []
        for line in lines[start_line:]:
            if 'Loop' in line:
                break
            data_lines.append(line)
    data = pd.read_csv(io.StringIO(''.join(data_lines)), delim_whitespace=True)
    # Create plots for each column (except 'CPU' and 'Step')
    columns_to_plot = [col for col in data.columns if col not in ['CPU', 'Step']]

    for col in columns_to_plot:
        plt.figure(figsize=(10, 6))
        plt.plot(data['Step']*dt/1000, data[col],color='black')
        plt.xlabel('Time [ps]')
        plt.ylabel(f'{col} [{get_units(col)}]')
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        plt.tick_params(direction='in', which='both', top=True, right=True)
        plt.tight_layout()
        plt.savefig(f'{col}.pdf', dpi=1000)
        plt.close()

def get_units(column):
    units = {
        'Temp': 'K',
        'PotEng': 'kcal/mol',
        'KinEng': 'kcal/mol',
        'TotEng': 'kcal/mol',
        'Density': 'g/cm³'
    }
    return units.get(column, '')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot LAMMPS log file properties')
    parser.add_argument('log_file', type=str, help='Path to the LAMMPS log file')
    parser.add_argument('-dt','--timestep', type=float, help='MD time step in fs')
    args = parser.parse_args()

    plot_lammps_log(args.log_file,args.timestep)
