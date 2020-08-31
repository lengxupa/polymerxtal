"""
Functions for manipulating LAMMPS files.
"""

import os


def read_data(data_file):
    """Read information from a LAMMPS data file

    Parameters
    ----------
    data_file : str
        The path of the data file to read

    Returns
    -------
    data_dir : dict
        Data from the LAMMPS data file as a dictionary.
    """

    # This function would be easier to unit test if it took in a string instead of a file.
    data_dir = {}
    box = {}
    masses = {}
    with open(data_file, 'r') as data_src:
        skip_line = 1
        feature = ''
        for line in data_src.readlines():
            if skip_line:
                skip_line = 0
                continue

            # ln will have 0 length if line is comment.
            ln = line.split('#')[0].split()
            if ln:
                if len(ln) > 1 and ln[0].isdigit() and (not ln[1].isdigit()) and (not feature):
                    # I think you are just trying to cast as a number here? Is safer to try casting as float.
                    data_dir[' '.join(ln[1:])] = float(ln[0])
                if len(ln) == 4 and ln[2][1:] == 'lo' and ln[3][1:] == 'hi':
                    data_dir[ln[2]] = float(ln[0])
                    data_dir[ln[3]] = float(ln[1])
                if not (ln[0][0].isdigit() or ln[0][0] == '-'):
                    feature = ' '.join(ln)
                    data_dir[feature] = {}
                if feature and (ln[0][0].isdigit() or ln[0][0] == '-'):
                    data_dir[feature][eval(ln[0])] = [eval(i) for i in ln[1:]]
    return data_dir


def write_lmp_ifile(file_location, datafile, potential_headfile, potentialfile):

    # Write a LAMMPS input file given a file location and datafile.
    # Where/what is in.lapps? Seems like a template.
    # I have removed in.lammps and changed to a direct write
    des = open('inr.lammps', 'w')
    des.write('newton          on\n')
    des.write('boundary        p p p\n')
    des.write('units           real\n')
    des.write('\n')
    des.write('include      ../%s\n' % potential_headfile)
    des.write('read_data       ../%s\n' % datafile)
    des.write('include      ../%s\n' % potentialfile)
    des.write('\n')
    des.write('compute         graincm all com\n')
    des.write('variable        M equal mass(all)\n')
    des.write('variable        maxcx equal bound(all,xmax)\n')
    des.write('variable        mincx equal bound(all,xmin)\n')
    des.write('variable        maxcy equal bound(all,ymax)\n')
    des.write('variable        mincy equal bound(all,ymin)\n')
    des.write('variable        maxcz equal bound(all,zmax)\n')
    des.write('variable        mincz equal bound(all,zmin)\n')
    des.write('\n')
    des.write(
        'fix             1 all ave/time 1 1 1 c_graincm[1] c_graincm[2] c_graincm[3] v_maxcx v_mincx v_maxcy v_mincy v_maxcz v_mincz v_M file tmp.out\n'
    )
    # Here is the tmp.out output after running LAMMPS
    des.write('run             0\n')
    des.write('\n')
    des.close()


def run_lammps(lmp_location):
    os.system('%s < inr.lammps' %lmp_location)
    os.system('cp log.lammps out')


def read_cell_sizes(data_file):
    # read_cell_sizes function provides a quick step to read the simulation box information
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
    # get_boundaries function gets the center of mass and the boundaris - the minimum and maximum x, y and z positions - of a molecule or molecules from LAMMPS out put tmp.out 
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
