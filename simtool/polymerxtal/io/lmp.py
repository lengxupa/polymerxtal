"""
Functions for manipulating LAMMPS files.
"""

import os
from subprocess import call, Popen, PIPE
from time import strftime

LAMMPS_EXEC = os.environ.get('LAMMPS_EXEC')
verbose = False


def check_lmps_exec():
    if LAMMPS_EXEC is None:
        print('you must set environment variable LAMMPS_EXEC')
        return False
    else:
        try:
            stdout, stderr = Popen([LAMMPS_EXEC, '-e', 'both', '-l', 'none'], stdin=PIPE, stdout=PIPE,
                                   stderr=PIPE).communicate()
            if verbose:
                print('using %s LAMMPS machine' % LAMMPS_EXEC)
            return True
        except OSError:
            print('LAMMPS is not configured properly for one reason or another')
            return False


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


def write_lmp_ifile(lmp_file_location, datafile, potential_headfile, potentialfile):

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


def call_lammps(lmp_input, np, nanohub, prefix='mpiexec', print_to_screen=False):
    """pysimm.lmps.call_lammps
    Wrapper to call LAMMPS using executable name defined in pysimm.lmps module.
    Args:
        simulation: :class:`~pysimm.lmps.Simulation` object reference
        np: number of threads to use
        nanohub: dictionary containing nanohub resource information default=None
        prefix: prefix for running LAMMPS (i.e. - mpiexec)
    Returns:
        None
    """

    log_name = 'out'

    #if nanohub:
    #    with open('temp.in', 'w') as f:
    #        f.write(simulation.input)
    #    if simulation.name:
    #        print('%s: sending %s simulation to computer cluster at nanoHUB' % (strftime('%H:%M:%S'), simulation.name))
    #    else:
    #        print('%s: sending simulation to computer cluster at nanoHUB' % strftime('%H:%M:%S'))
    #    sys.stdout.flush()
    #    cmd = ('submit -n %s -w %s -i temp.lmps -i temp.in '
    #           'lammps-09Dec14-parallel -e both -l none -i temp.in'
    #           % (nanohub.get('cores'), nanohub.get('walltime')))
    #    cmd = shlex.split(cmd)
    #    exit_status, stdo, stde = RapptureExec(cmd)
    #else:
    #    if simulation.name:
    #        print('%s: starting %s LAMMPS simulation'
    #              % (strftime('%H:%M:%S'), simulation.name))
    #    else:
    #        print('%s: starting LAMMPS simulation'
    #              % strftime('%H:%M:%S'))
    print('%s: starting LAMMPS simulation' % strftime('%H:%M:%S'))
    if np:
        p = Popen([prefix, '-np', str(np), LAMMPS_EXEC, '-e', 'both'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        p = Popen([LAMMPS_EXEC, '-e', 'both'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    #    simulation.write_input()
    #    if simulation.debug:
    #        print(simulation.input)
    #        warning_print('debug setting involves streaming output from LAMMPS process and can degrade performance')
    #        warning_print('only use debug for debugging purposes, use print_to_screen to collect stdout after process finishes')
    #        p.stdin.write(simulation.input)
    #        q = Queue()
    #        t = Thread(target=enqueue_output, args=(p.stdout, q))
    #        t.daemon = True
    #        t.start()
    #
    #        while t.isAlive() or not q.empty():
    #            try:
    #                line = q.get_nowait()
    #            except Empty:
    #                pass
    #            else:
    #                if simulation.debug:
    #                    sys.stdout.write(line)
    #                    sys.stdout.flush()
    #    else:
    simulation_input = ""
    with open(lmp_input, 'r') as input_text:
        simulation_input = input_text.read()
    stdo, stde = p.communicate(simulation_input.encode('utf-8'))
    if print_to_screen:
        print(stdo)
        print(stde)
    if stde:
        return 1
    else:
        return 0

    #simulation.system.read_lammps_dump('pysimm.dump.tmp')

    #try:
    #    os.remove('temp.lmps')
    #except OSError as e:
    #    print(str(e))

    #if os.path.isfile('pysimm.qeq.tmp'):
    #    os.remove('pysimm.qeq.tmp')

    #try:
    #    os.remove('pysimm.dump.tmp')
    #    if simulation.name:
    #        print('%s: %s simulation using LAMMPS successful'
    #              % (strftime('%H:%M:%S'), simulation.name))
    #    else:
    #        print('%s: simulation using LAMMPS successful'
    #              % (strftime('%H:%M:%S')))
    #except OSError as e:
    #    if simulation.name:
    #        raise PysimmError('%s simulation using LAMMPS UNsuccessful' % simulation.name)
    #    else:
    #        raise PysimmError('simulation using LAMMPS UNsuccessful')


def run_lammps(lmp_input, np=None, nanohub=None, save_input=True, prefix='mpiexec'):
    """pysimm.lmps.Simulation.run
    Begin LAMMPS simulation.
    Args:
        np: number of threads to use (serial by default) default=None
        nanohub: dictionary containing nanohub resource information default=None
        init: True to write initialization part of LAMMPS input script (set to False if using complete custom input)
        save_input: True to save input as pysimm.sim.in
        prefix: prefix for running LAMMPS (i.e. - mpiexec)
    """
    try:
        return_code = call_lammps(lmp_input, np, nanohub, prefix=prefix)
    except OSError:
        raise Exception('There was a problem calling LAMMPS with {}'.format(prefix))
    except IOError:
        if check_lmps_exec():
            raise Exception(
                'There was a problem running LAMMPS. The process started but did not finish successfully. Check the log file, or rerun the simulation with debug=True to debug issue from LAMMPS output'
            )
        else:
            raise Exception(
                'There was a problem running LAMMPS. LAMMPS is not configured properly. Make sure the LAMMPS_EXEC environment variable is set to the correct LAMMPS executable path. The current path is set to:\n\n{}'
                .format(LAMMPS_EXEC))

    #os.system('%s < %s' % (lmp_location, lmp_file_location))
    #os.system('cp log.lammps out')

    #return_code = os.system('lmp_serial < ./bonds/bondcreate.in > out')
    return return_code

    def run(self, np=None, nanohub=None, save_input=True, prefix='mpiexec'):
        """pysimm.lmps.Simulation.run
        Begin LAMMPS simulation.
        Args:
            np: number of threads to use (serial by default) default=None
            nanohub: dictionary containing nanohub resource information default=None
            init: True to write initialization part of LAMMPS input script (set to False if using complete custom input)
            save_input: True to save input as pysimm.sim.in
            prefix: prefix for running LAMMPS (i.e. - mpiexec)
        """
        if isinstance(save_input, str):
            with open(save_input, 'w') as f:
                f.write(self.input)
        elif save_input is True:
            with open('pysimm.sim.in', 'w') as f:
                f.write(self.input)
        try:
            call_lammps(self, np, nanohub, prefix=prefix)
        except OSError:
            raise PysimmError('There was a problem calling LAMMPS with {}'.format(prefix))
        except IOError:
            if check_lmps_exec():
                raise PysimmError(
                    'There was a problem running LAMMPS. The process started but did not finish successfully. Check the log file, or rerun the simulation with debug=True to debug issue from LAMMPS output'
                )
            else:
                raise PysimmError(
                    'There was a problem running LAMMPS. LAMMPS is not configured properly. Make sure the LAMMPS_EXEC environment variable is set to the correct LAMMPS executable path. The current path is set to:\n\n{}'
                    .format(LAMMPS_EXEC))


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
