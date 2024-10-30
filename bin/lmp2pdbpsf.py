#!/usr/bin/env python
# $Id$

__doc__ = """
Reads a LAMMPS data file and outputs a pdb/psf combo
"""

import string, getopt, sys, os

try:
    import numpy as N
except ImportError:
    try:
        import Numeric as N
    except ImportError:
        print "Error: You must have either Numerical Python or numpy installed"
        sys.exit(-1)

#------------------------------------------------------------
# Classes
#------------------------------------------------------------

class Atom(object):
    __slots__ = ("index", "name", "type", "chainid", "charge", "mass", "chain_index", "_positions")
    def __init__(self, index, name, type, chain_id, charge=0, mass=1):
        self.index = index
        self.name = repr(type)
        self.type = type
        self.chainid = chain_id
        self.charge = charge
        self.mass = mass
    def __repr__(self):
        return "<Atom "+repr(self.index+1)+ ": name " + repr(self.type) +" of chain "+repr(self.chainid)+">"
    def __cmp__(self, other):
        return cmp(self.index, other.index)
    def __eq__(self, other):
        return self.index == other.index
    def __hash__(self):
        return hash(self.index)
    def __getattr__(self, attr):
        if attr == 'pos':
            return self._positions[self.index]
        else: super(Atom, self).__getattribute__(attr)

class LAMMPSParseError(Exception):
    pass

header_keywords= ["atoms","bonds","angles","dihedrals","impropers","atom types","bond types","angle types","dihedral types","improper types","xlo xhi","ylo yhi","zlo zhi"]
connections = dict([["Bonds",("bonds", 3)],["Angles",("angles", 3)],
                    ["Dihedrals",("dihedrals", 4)],["Impropers",("impropers", 2)]])
coeff = dict([["Masses",("atom types", 1)], ["Velocities",("atoms", 3)],
             ["Pair Coeffs",("atom types", 4)],
             ["Bond Coeffs",("bond types", 2)],["Angle Coeffs",("angle types", 4)],
             ["Dihedral Coeffs",("dihedral types", 3)], ["Improper Coeffs",("improper types", 2)]])

def conv_float(l):
    """
    Function to be passed into map or a list comprehension. If the argument is a float it is converted,
    otherwise the original string is passed back
    """
    try:
        n = float(l)
    except ValueError:
        n = l
    return n

class LAMMPSData(object):
    def __init__(self, filename=None):
        self.names = {}
        self.headers = {}
        self.sections = {}
        if filename == None:
            self.title = "LAMMPS data file"
        else:
            # Open and check validity
            file = open(filename, 'r')
            file_iter = file.xreadlines()
            self.title = file_iter.next()
            # Parse headers
            headers = self.headers
            for l in file_iter:
                line = l.strip() 
                if len(line) == 0: continue
                found = False
                for keyword in header_keywords:
                    if line.find(keyword) >= 0:
                        found = True
                        values = line.split()
                        if keyword in ("xlo xhi", "ylo yhi", "zlo zhi"):
                            headers[keyword] = (float(values[0]), float(values[1]))
                        else:
                            headers[keyword] = int(values[0])
                if found == False: break
            file.close()
            # Parse sections
            # XXX This is a crappy way to do it
            file = open(filename, 'r')
            file_iter = file.xreadlines()
            # Create coordinate array
            positions = N.zeros((headers['atoms'], 3), N.Float64)
            sections = self.sections
            for l in file_iter:
                line = l.strip()
                if len(line) == 0: continue
                if coeff.has_key(line):
                    h, numcoeff = coeff[line]
                    # skip line
                    file_iter.next()
                    data = []
                    for i in xrange(headers[h]):
                        fields = file_iter.next().strip().split()
                        data.append(tuple(map(conv_float, fields[1:])))
                    sections[line] = data
                elif connections.has_key(line):
                    h, numfields = connections[line]
                    # skip line
                    file_iter.next()
                    data = []
                    for i in range(headers[h]):
                        fields = file_iter.next().strip().split()
                        data.append(tuple(map(int, fields[1:])))
                    sections[line] = data
                elif line == "Atoms":
                    file_iter.next()
                    data = []
                    for i in xrange(headers["atoms"]):
                        fields = file_iter.next().strip().split()
                        a = Atom(index=int(fields[0])-1, name=fields[2], type=int(fields[2])-1, chain_id=int(fields[1]), charge=float(fields[3]))
                        a._positions = positions
                        data.append(a)
                        positions[i] = N.array([float(fields[4]), float(fields[5]), float(fields[6])])
                    sections[line] = data
            self.positions = positions
            file.close()
    def writePSF(self, filename):
        import string
        psf_atom_format = "   %5d %4s %4d %4s %-4s %-4s %10.6f      %7.4f            %1d\n"
        file = open(filename, 'w')
        file.write("PSF\n")
        file.write(string.rjust(str(len(self.sections["Atoms"])), 8) + ' !NATOM\n')
        for i, atom in enumerate(self.sections["Atoms"]):
            line = [i+1, 'TEMP', atom.chainid, 'XXXX', str(atom.type+1), str(atom.type+1), 0.0, 0.0, 0]
            file.write(psf_atom_format%tuple(line))
        file.write("\n")
        num_bonds = len(self.sections["Bonds"])
        bond_list = self.sections["Bonds"]
        file.write(string.rjust(str(num_bonds), 8) + ' !NBOND\n')
        for index in range(0, num_bonds, 4):
            try:
                bonds = bond_list[index:index+4]
            except IndexError:
                bonds = bond_list[index:-1]
            bond_line = map(lambda bond: string.rjust(str(bond[1]), 8)+string.rjust(str(bond[2]), 8), bonds)
            file.write(''.join(bond_line)+'\n')
        file.close()
    def writePDB(self, filename):
        import string
        atom_format = "%6s%5d %4s %4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
        p = self.positions
        file = open(filename, 'w')
        for i, atom in enumerate(self.sections["Atoms"]):
            line = ["ATOM  ", i+1, str(atom.chainid), 'XXXX', atom.type+1, p[i,0], p[i,1], p[i,2], 0.0, 0.0, str(atom.type)]
            file.write(atom_format%tuple(line))


#------------------------------------------------------------
# functions
#------------------------------------------------------------

def usage(rc):
    """
       Print usage to stdout.

       Parameters
            rc:  Exit status (int)
    """
    
    str = """usage: %s [OPTIONS] -d in.lammps -o out

Read a LAMMPS data file and write out a psf/pdb combo from the information in the file. 
This is useful if you wish to view the system in VMD or other molecular viewer. You can 
output just the psf, which is usefull to view LAMMPS trajectory files (dump or dcd)

-h,--help             this help
-d,--lammps=FILE      lammps input file
-o,--prefix=FILE      prefix for pdb/psf output
-p,--psf              only output PSF file
-s,--pdb              only output PDB file
""" % sys.argv[0]
    sys.stderr.write(str)
    sys.exit(rc)

def mainCommand():
    """
        Main driver for running program from the command line.
    """
    shortOptlist = "hd:o:p:s"
    longOptlist = ["help","lammps=","prefix=","psf","pdb"]
    if len(sys.argv) == 1: usage(0)

    try: opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)

    infile = None
    outfile_prefix = None
    only_pdb = False
    only_psf = False
    for o,a in opts:
        if o in ("-h","--help"):
            usage(0)
        elif o in ("-d","--lammps"):
            infile = a
        elif o in ("-o","--prefix"):
            outfile_prefix = a
        elif o in ("-p","--psf"):
            only_psf = True
        elif o in ("-s","--pdb"):
            only_pdb = True
    if only_pdb == True and only_psf == True:
        sys.stderr.write("Option error: Only one exclusive option can be specified\n")
        sys.exit(-1)

    if infile == None:
        sys.stderr.write("Error: must specify lammps data file\n")
        sys.exit(-1)
    if outfile_prefix == None:
        outfile_prefix = os.path.splitext(os.path.split(infile)[1])[0]
    d = LAMMPSData(infile)
    if only_pdb and only_psf: only_pdb = only_psf = False
    if not only_pdb: d.writePSF(outfile_prefix+".psf")
    if not only_psf: d.writePDB(outfile_prefix+".pdb")

if __name__ == "__main__":
    mainCommand()
