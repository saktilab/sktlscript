#!/usr/bin/env python3

"""
Calculate RMSD between two XYZ files


"""

import numpy as np
import re


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P = rotate(P, Q)
    return rmsd(P, Q)


def rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def write_coordinates(atoms, V):
    """
    Print coordinates V
    """
    N, D = V.shape

    print (str(N))
    print()

    for i in range(N):
        line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(atoms[i], V[i, 0], V[i, 1], V[i, 2])
        print (line)

def get_coordinates(filename, fmt, ignore_hydrogens):
    """
    Get coordinates from filename.

    """
    if fmt == "xyz":
        return get_coordinates_xyz(filename, ignore_hydrogens)
    elif fmt == "pdb":
        return get_coordinates_pdb(filename, ignore_hydrogens)
    exit("Could not recognize file format: {:s}".format(fmt))


def get_coordinates_pdb(filename, ignore_hydrogens):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.

    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    V = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the atomtype
                try:
                    atom = tokens[2][0]
                    if ignore_hydrogens and atom == "H":
                        continue
                    elif atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                except:
                        exit("Error parsing atomtype for the following line: \n%s" % line)

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n%s" % line)
                # Try to read the coordinates
                try:
                    V.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        V.append(np.asarray([x,y,z],dtype=float))
                    except:
                        exit("Error parsing input for the following line: \n%s" % line)


    V = np.asarray(V)
    return atoms, V


def get_coordinates_xyz(filename, ignore_hydrogens):
    """
    Get coordinates from a filename.xyz and return a vectorset with all the
    coordinates.

    This function has been written to parse XYZ files, but can easily be
    written to parse others.

    """
    f = open(filename, 'r')
    V = []
    atoms = []
    n_atoms = 0
    lines_read = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(next(f))
    except ValueError:
        exit("Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    next(f)

    # Use the number of atoms to not read beyond the end of a file
    for line in f:

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        # ignore hydrogens
        if ignore_hydrogens and atom.lower() == "h":
            continue

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) == 3:
            V.append(np.array(numbers))
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

        lines_read += 1

    f.close()
    V = np.array(V)
    return atoms, V


if __name__ == "__main__":
    import argparse

    description = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ or PDB format.
The order of the atoms *must* be the same for both structures.
"""

    epilog = """
The script will return three RMSD values:
1) Normal: The RMSD calculated the straight-forward way.
2) Kabsch: The RMSD after the two coordinate sets are translated and rotated onto each other.
"""

    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('structure_a', metavar='structure_a.xyz', type=str)
    parser.add_argument('structure_b', metavar='structure_b.xyz', type=str)
    parser.add_argument('-o', '--output', action='store_true', help='print out structure A, centered and rotated unto structure B\'s coordinates in XYZ format')
    parser.add_argument('-n', '--no-hydrogen', action='store_true', help='ignore hydrogens when calculating RMSD')
    parser.add_argument('-f', '--format', action='store', help='Format of input files. Supports xyz or pdb.', default="xyz")

    args = parser.parse_args()

    atomsP, P = get_coordinates(args.structure_a, args.format, args.no_hydrogen)
    atomsQ, Q = get_coordinates(args.structure_b, args.format, args.no_hydrogen)

    # Calculate 'dumb' RMSD
    normal_rmsd = rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    if args.output:
        V = rotate(P, Q)
        V += Qc
        write_coordinates(atomsP, V)
        quit()

    print ("Normal RMSD:", normal_rmsd)
    print ("Kabsch RMSD:", kabsch_rmsd(P, Q))
