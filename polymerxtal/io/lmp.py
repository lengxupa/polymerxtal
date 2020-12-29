"""
Functions for manipulating LAMMPS files.
"""

import os
from subprocess import call, Popen, PIPE
from time import strftime

from .hublib import check_nanohub

LAMMPS_EXEC = os.environ.get("LAMMPS_EXEC")
if LAMMPS_EXEC is None:
    LOADEDMODULES = os.environ.get("LOADEDMODULES")
    if LOADEDMODULES and ("lammps/31Mar17" in LOADEDMODULES):
        LAMMPS_EXEC = "lmp"
verbose = False


def check_lmps_exec():
    """Check lammps executable exists"""
    if LAMMPS_EXEC is None:
        print("you must set environment variable LAMMPS_EXEC")
        return False
    else:
        try:
            stdout, stderr = Popen(
                [LAMMPS_EXEC, "-e", "both", "-l", "none"],
                stdin=PIPE,
                stdout=PIPE,
                stderr=PIPE,
            ).communicate()
            if verbose:
                print("using %s LAMMPS machine" % LAMMPS_EXEC)
            return True
        except OSError:
            print("LAMMPS is not configured properly for one reason or another")
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
    with open(data_file, "r") as data_src:
        skip_line = 1
        feature = ""
        for line in data_src.readlines():
            if skip_line:
                skip_line = 0
                continue

            # ln will have 0 length if line is comment.
            ln = line.split("#")[0].split()
            if ln:
                if (
                    len(ln) > 1
                    and ln[0].isdigit()
                    and (not ln[1].isdigit())
                    and (not feature)
                ):
                    # I think you are just trying to cast as a number here? Is safer to try casting as float.
                    data_dir[" ".join(ln[1:])] = float(ln[0])
                if len(ln) == 4 and ln[2][1:] == "lo" and ln[3][1:] == "hi":
                    data_dir[ln[2]] = float(ln[0])
                    data_dir[ln[3]] = float(ln[1])
                if not (ln[0][0].isdigit() or ln[0][0] == "-"):
                    feature = " ".join(ln)
                    data_dir[feature] = {}
                if feature and (ln[0][0].isdigit() or ln[0][0] == "-"):
                    data_dir[feature][eval(ln[0])] = [eval(i) for i in ln[1:]]
    return data_dir


def write_lmp_ifile(
    lmp_file_location="LAMMPSinputfile.txt",
    datafile="",
    potential_headfile="",
    potentialfile="",
):
    """Write a LAMMPS input file given a file location and datafile.

    Parameters
    ----------
    lmp_file_location : str (optional)
        The path of LAMMPS input file to output.
    datafile : str (optional)
        The path of the LAMMPS data file to read.
    potential_headfile : str (optional)
        The path of part of the LAMMPS input file including potential style definitions.
    potentialfile : str (optional)
        The path of part of the LAMMPS input file including potential coefficients.

    Returns
    -------
    None
    """

    des = open(lmp_file_location, "w")
    # des.write("newton          on\n")
    des.write("boundary        p p p\n")
    des.write("units           real\n")
    des.write("\n")
    des.write("atom_style        full\n")
    des.write("special_bonds     lj/coul 0.0 0.0 1.0 dihedral yes\n")
    des.write("dielectric        1.0\n")
    des.write("pair_style        lj/cut  12.0\n")
    des.write("bond_style        harmonic\n")
    des.write("angle_style       harmonic\n")
    des.write("dihedral_style    harmonic\n")
    des.write("improper_style    harmonic\n")
    if potential_headfile:
        des.write("include      %s\n" % potential_headfile)
    if datafile:
        des.write("read_data       %s\n" % datafile)
    des.write("neighbor          0.3 bin\n")
    des.write(
        "thermo_style      custom step etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press pxx pyy pzz pxy pxz pyz lx ly lz vol density\n"
    )
    des.write("thermo            10\n")
    des.write("thermo_modify     flush yes\n")
    des.write("\n")
    des.write("pair_style        buck/coul/long  12.0 12.0\n")
    des.write("kspace_style      pppm 1e-4\n")
    if potentialfile:
        des.write("include      %s\n" % potentialfile)
    des.write("\n")
    des.write("fix 1 all nve\n")
    des.write("run               0\n")
    des.write("\n")
    des.write("# Minimization parameters\n")
    des.write("minimize          1.0e-9 1.0e-9 5000 100000\n")
    des.write("# Dump minimized system\n")
    des.write("dump              1 all atom 1 min.dump\n")
    des.write("dump_modify       1 image yes scale no\n")
    des.write("fix 1 all nve\n")
    des.write("run               0\n")
    des.write("unfix 1\n")
    des.write("\n")
    des.write("undump            1\n")
    des.write("write_data min.data\n")
    des.write("\n")
    des.write("# MD parameters\n")
    des.write("run_style         respa 3 2 2 bond 1 pair 2 kspace 3\n")
    des.write("reset_timestep    0\n")
    des.write("timestep          4\n")
    des.write("dump              1 all atom 100 md.dump\n")
    des.write("dump_modify       1 image yes scale no\n")
    des.write("fix 1 all nvt temp 300.0 300.0 400.0\n")
    des.write("run 1000\n")
    des.write("write_restart restart.lammps\n")
    des.write("write_data md.data\n")
    des.write("\n")
    # des.write("compute         graincm all com\n")
    # des.write("variable        M equal mass(all)\n")
    # des.write("variable        maxcx equal bound(all,xmax)\n")
    # des.write("variable        mincx equal bound(all,xmin)\n")
    # des.write("variable        maxcy equal bound(all,ymax)\n")
    # des.write("variable        mincy equal bound(all,ymin)\n")
    # des.write("variable        maxcz equal bound(all,zmax)\n")
    # des.write("variable        mincz equal bound(all,zmin)\n")
    # des.write("\n")
    # des.write(
    #    "fix             1 all ave/time 1 1 1 c_graincm[1] c_graincm[2] c_graincm[3] v_maxcx v_mincx v_maxcy v_mincy v_maxcz v_mincz v_M file tmp.out\n"
    # )
    # Here is the tmp.out output after running LAMMPS
    # des.write("run             0\n")
    # des.write("\n")
    des.close()


def call_lammps(lmp_input, np, nanohub=None, prefix="mpiexec", print_to_screen=False):
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

    log_name = "out"

    if nanohub:
        # with open('temp.in', 'w') as f:
        #    f.write(simulation.input)
        # if simulation.name:
        #    print('%s: sending %s simulation to computer cluster at nanoHUB' % (strftime('%H:%M:%S'), simulation.name))
        # else:
        #    print('%s: sending simulation to computer cluster at nanoHUB' % strftime('%H:%M:%S'))
        # sys.stdout.flush()
        # cmd = ('submit -n %s -w %s -i temp.lmps -i temp.in '
        #       'lammps-09Dec14-parallel -e both -l none -i temp.in'
        #       % (nanohub.get('cores'), nanohub.get('walltime')))
        # cmd = shlex.split(cmd)
        # exit_status, stdo, stde = RapptureExec(cmd)
        return_code = os.system("lmp_serial < %s > out" % lmp_input)
        return return_code
    else:
        #    if simulation.name:
        #        print('%s: starting %s LAMMPS simulation'
        #              % (strftime('%H:%M:%S'), simulation.name))
        #    else:
        #        print('%s: starting LAMMPS simulation'
        #              % strftime('%H:%M:%S'))
        print("%s: starting LAMMPS simulation" % strftime("%H:%M:%S"))
        if np:
            p = Popen(
                [prefix, "-np", str(np), LAMMPS_EXEC, "-e", "both"],
                stdin=PIPE,
                stdout=PIPE,
                stderr=PIPE,
            )
        else:
            p = Popen([LAMMPS_EXEC, "-e", "both"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        # simulation.write_input()
        # if simulation.debug:
        #    print(simulation.input)
        #    warning_print('debug setting involves streaming output from LAMMPS process and can degrade performance')
        #    warning_print('only use debug for debugging purposes, use print_to_screen to collect stdout after process finishes')
        #    p.stdin.write(simulation.input)
        #    q = Queue()
        #    t = Thread(target=enqueue_output, args=(p.stdout, q))
        #    t.daemon = True
        #    t.start()

        #    while t.isAlive() or not q.empty():
        #        try:
        #            line = q.get_nowait()
        #        except Empty:
        #            pass
        #        else:
        #            if simulation.debug:
        #                sys.stdout.write(line)
        #                sys.stdout.flush()
        # else:
        simulation_input = ""
        with open(lmp_input, "r") as input_text:
            simulation_input = input_text.read()
        stdo, stde = p.communicate(simulation_input.encode("utf-8"))
        if print_to_screen:
            print(stdo)
            print(stde)
        if stde:
            return 1
        else:
            return 0

    # simulation.system.read_lammps_dump('pysimm.dump.tmp')

    # try:
    #    os.remove('temp.lmps')
    # except OSError as e:
    #    print(str(e))

    # if os.path.isfile('pysimm.qeq.tmp'):
    #    os.remove('pysimm.qeq.tmp')

    # try:
    #    os.remove('pysimm.dump.tmp')
    #    if simulation.name:
    #        print('%s: %s simulation using LAMMPS successful'
    #              % (strftime('%H:%M:%S'), simulation.name))
    #    else:
    #        print('%s: simulation using LAMMPS successful'
    #              % (strftime('%H:%M:%S')))
    # except OSError as e:
    #    if simulation.name:
    #        raise PysimmError('%s simulation using LAMMPS UNsuccessful' % simulation.name)
    #    else:
    #        raise PysimmError('simulation using LAMMPS UNsuccessful')


def run_lammps(lmp_input, np=None, nanohub=None, save_input=True, prefix="mpiexec"):
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
        return_code = call_lammps(lmp_input, np, nanohub=check_nanohub(), prefix=prefix)
    except OSError:
        raise Exception("There was a problem calling LAMMPS with {}".format(prefix))
    except:
        if check_lmps_exec():
            raise Exception(
                "There was a problem running LAMMPS. The process started but did not finish successfully. Check the log file, or rerun the simulation with debug=True to debug issue from LAMMPS output"
            )
        else:
            raise Exception(
                "There was a problem running LAMMPS. LAMMPS is not configured properly. Make sure the LAMMPS_EXEC environment variable is set to the correct LAMMPS executable path. The current path is set to:\n\n{}".format(
                    LAMMPS_EXEC
                )
            )

    # os.system('%s < %s' % (lmp_location, lmp_file_location))
    # os.system('cp log.lammps out')

    # return_code = os.system('lmp_serial < ./bonds/bondcreate.in > out')
    return return_code


def read_cell_sizes(data_file):
    # read_cell_sizes function provides a quick step to read the simulation box information
    cell_sizes = {}
    delta_cell = {}
    with open(data_file) as f:
        for line in f:
            if line.endswith("xhi\n"):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                ds = float(hi) - float(lo)
                cell_sizes["x"] = [lo, hi, ds]
                delta_cell["x"] = ds
            elif line.endswith("yhi\n"):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                ds = float(hi) - float(lo)
                cell_sizes["y"] = [lo, hi, ds]
                delta_cell["y"] = ds
            elif line.endswith("zhi\n"):
                values = line.rstrip().split()
                lo = values[0]
                hi = values[1]
                ds = float(hi) - float(lo)
                cell_sizes["z"] = [lo, hi, ds]
                delta_cell["z"] = ds

                break

    return cell_sizes, delta_cell


def get_boundaries(direction):
    # get_boundaries function gets the center of mass and the boundaris - the minimum and maximum x, y and z positions - of a molecule or molecules from LAMMPS out put tmp.out
    boundaries = {}
    com = {}
    with open("tmp.out", "r") as infile:
        r = infile.readline()
        r = infile.readline()
        for line in infile:
            print(line)
            line = [float(x) for x in line.split()]
    mass = line[10]
    com["x"] = line[1]
    com["y"] = line[2]
    com["z"] = line[3]
    boundaries["x"] = [line[4], line[5]]
    boundaries["y"] = [line[6], line[7]]
    boundaries["z"] = [line[8], line[9]]
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
        # for v in values:
        #    distance.append(abs(v-com[key]))

    radius = math.ceil(max(distance)) + 5

    #    print('Coordinates cylinder: ', coord, 'Radius distance :', radius)
    return mass, com, boundaries, coord, radius


def write_data(data_dir, ofile, velocities=True, atom_only=False):
    """Write a LAMMPS data file

    Parameters
    ----------
    dir_data : dict
        Dictionary storing data to form LAMMPS data file.
    ofile : str
        The path of LAMMPS data file to output.
    velocities : bool (optional)
        Decide ouput velocities information or not.
    atom_only : bool (optional)
        Decide output atoms information only or not.

    Returns
    -------
    data_dir : dict
        Data from the LAMMPS data file as a dictionary.
    """
    des = open(ofile, "w")
    des.write("LAMMPS data file via Tongtong\n")
    des.write("\n")
    if atom_only:
        head_list = ["atom"]
        main_list = ["Masses", "Atoms", "Velocities"]
    else:
        head_list = ["atom", "bond", "angle", "dihedral", "improper"]
        main_list = [
            "Masses",
            "Pair Coeffs",
            "Bond Coeffs",
            "Angle Coeffs",
            "Dihedral Coeffs",
            "Improper Coeffs",
            "Atoms",
            "Velocities",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
        ]
    for head in head_list:
        if (head + "s") in data_dir:
            des.write("%d %s\n" % (data_dir[head + "s"], (head + "s")))
            des.write("%d %s\n" % (data_dir[head + " types"], (head + " types")))
    des.write("\n")
    for xyz in ["x", "y", "z"]:
        des.write(
            "%f %f %s %s\n"
            % (data_dir[xyz + "lo"], data_dir[xyz + "hi"], (xyz + "lo"), (xyz + "hi"))
        )
    des.write("\n")
    if not velocities:
        main_list.remove("Velocities")
    for key in main_list:
        if key in data_dir and len(data_dir[key]):
            des.write(key + (" # full\n" if key == "Atoms" else "\n"))
            des.write("\n")
            for i in data_dir[key]:
                des.write(
                    str(i) + " " + " ".join(str(j) for j in data_dir[key][i]) + "\n"
                )
            des.write("\n")
    des.close()
