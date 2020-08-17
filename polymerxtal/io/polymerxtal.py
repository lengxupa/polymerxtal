"""
Functions for manipulating PolymerXtal files.
"""


def read_input(infile):
    ifile = open(infile)
    args = {}
    for line in ifile.readlines():
        valid_line = line.split('#')[0]
        ln = valid_line.split()
        if len(ln) < 1:
            continue
        args[ln[0]] = ln[1:]
    ifile.close()
    return args


def check_args(args):
    # Check data file input
    if 'datafile' not in args:
        print(
            'Please specify HE grain LAMMPS datafile, (please refer to data.1.in and data.6.in in examples directory for format)'
        )
        return False

    # Check potential header file input
    if 'potential_headfile' not in args:
        print(
            'Please specify HE grain LAMMPS potential head file, (please refer to potential_head.mod in examples directory for format)'
        )
        return False

    # Check potential file input
    if 'potentialfile' not in args:
        print(
            'Please specify HE grain LAMMPS potential file, (please refer to potential.mod in examples directory for format)'
        )
        return False
