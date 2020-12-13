"""
Functions for manipulating pdb files.
"""

import numpy as np


def open_pdb(file_location):
    """This function reads in a pdb file and returns the atom names and coordinates.

    The pdb file must specify the atom elements in the last column, and follow
    the conventions outlined in the PDB format specification.
    
    """

    with open(file_location) as f:
        data = f.readlines()

    coordinates = []
    symbols = []
    for line in data:
        if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
            symbols.append(line[76:79].strip())
            atom_coords = [float(x) for x in line[30:55].split()]
            coordinates.append(atom_coords)

    coords = np.array(coordinates)
    symbols = np.array(symbols)

    return symbols, coords


def write_pdb(file_location, symbols, coordinates):
    i_atom = 1

    RESIDUE_NAME = "UNK"
    CHAIN_ID = "A"
    RESIDUE_SEQ = 1

    # Write a pdb file given a file location, symbols, and coordinates.
    num_atoms = len(symbols)

    with open(file_location, 'w+') as f:
        #f.write('{}\n'.format(num_atoms))
        #f.write('XYZ file\n')

        for i in range(num_atoms):
            pos = coordinates[i]
            #        1     7   13  18     23     31   39   47
            f.write("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f\n" %
                    (i_atom, symbols[i], RESIDUE_NAME, CHAIN_ID, RESIDUE_SEQ, pos[0], pos[1], pos[2]))
            i_atom += 1
            if i_atom > 99999:
                i_atom = 1
        #        1     7
        f.write("TER   %5d\n" % i_atom)
        i_atom += 1
