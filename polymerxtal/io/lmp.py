"""
Functions for manipulating LAMMPS files.
"""

import os


def write_lmp_ifile(file_location, datafile):

    # Write a LAMMPS input file given a file location and datafile.
    with open('in.lammps', 'r') as file:
        filedata = file.read()

    filedata = filedata.replace('datafile', datafile)

    with open(file_location, 'w+') as file:
        file.write(filedata)


def run_lammps():
    os.system('./lmp_mpi < inr.lammps')
    os.system('cp log.lammps out')


def read_cell_sizes(data_file):
    cell_sizes = {}
    delta_cell = {}
    with open(data_file) as f:
        for line in f:
            if line.endswith('xhi\n'):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                ds = float(hi) - float(lo)
                cell_sizes['x'] = [lo, hi, ds]
                delta_cell['x'] = ds
            elif line.endswith('yhi\n'):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                ds = float(hi) - float(lo)
                cell_sizes['y'] = [lo, hi, ds]
                delta_cell['y'] = ds
            elif line.endswith('zhi\n'):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                ds = float(hi) - float(lo)
                cell_sizes['z'] = [lo, hi, ds]
                delta_cell['z'] = ds

                break

    return cell_sizes, delta_cell


def get_boundaries(direction):
    boundaries = {}
    com = {}
    with open('tmp.out', 'r') as infile:
        r = infile.readline()
        r = infile.readline()
        for line in infile:
            print(line)
            line = [float(x) for x in line.split()]
    mass = line[10]
    com['x'] = line[1]
    com['y'] = line[2]
    com['z'] = line[3]
    boundaries['x'] = [line[4], line[5]]
    boundaries['y'] = [line[6], line[7]]
    boundaries['z'] = [line[8], line[9]]
    # Get max radius
    distance = []
    coord = {}
    polymod = {}
    for key, values in sorted(boundaries.items()):
        if key == direction:
            continue
        res = sum(values) / 2.0
        coord[key] = res
        diff = abs(values[0] - res)
        distance.append(diff)
        #for v in values:
        #    distance.append(abs(v-com[key]))

    radius = math.ceil(max(distance)) + 5

    #    print('Coordinates cylinder: ', coord, 'Radius distance :', radius)
    return mass, com, boundaries, coord, radius
