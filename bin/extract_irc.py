#!/usr/bin/env python3

def extract_irc_coordinates(gaussian_output_file):
    coordinates = []
    with open(gaussian_output_file, 'r') as file:
        lines = file.readlines()

    # Find the line number where IRC coordinates start
    irc_start_line = None
    for i, line in enumerate(lines):
        if 'Coordinates' in line:
            irc_start_line = i + 5
            break

    if irc_start_line is None:
        print("No IRC coordinates found in the Gaussian output file.")
        return []

    # Extract the IRC coordinates
    for line in lines[irc_start_line:]:
        line = line.strip()
        if not line:
            break
        split_line = line.split()
        coordinates.append([float(x) for x in split_line[1:]])

    return coordinates

# Example usage
gaussian_output_file = 'inp.log'  # Replace with your Gaussian output file name
irc_coordinates = extract_irc_coordinates(gaussian_output_file)
print(irc_coordinates)

