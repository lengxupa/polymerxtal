"""
functions.py
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal
structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer,
tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given
atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers
where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline
polymer simulations.

Handles the primary functions
"""

import os, sys, os.path  # , math  # noqa: E401

# import shutil
import numpy as np


def calc_polymer(mass, nm, polymer_type, p_fra, custom=0):
    p = polymer_type.split("-")[0]
    wunits = {
        "PS": 104.1,
        "PE": 28.0516,
        "PMMA": 100.117,
        "PP": 42.0804,
        "POM": 30.0262,
        "PTFE": 100.0156,
        "PVC": 62.4987,
    }
    if p == "custom":
        wunit = custom
    else:
        wunit = wunits[p]
    m_polymer = (mass * p_fra) / (1 - p_fra)
    dp = int(m_polymer / wunit)
    nc = dp / nm
    print("Degree of polymerization: ", dp, "Chains: ", nc)
    return nc


def parameters(periodic, cell, coord, boundaries, radius):
    shift = {}

    if periodic == "x":
        lx = cell["x"][2]
        ly = (boundaries["y"][0] + radius) - (boundaries["y"][1] - radius)
        lz = (boundaries["z"][0] + radius) - (boundaries["z"][1] - radius)
        shift["y"] = boundaries["y"][1] - radius
        shift["z"] = boundaries["z"][1] - radius
        dx = 1
        dy = 0
        dz = 0
        # cx = -lx
        # cy = coord["y"] - shift["y"]
        # cz = coord["z"] - shift["z"]
        lc = lx
        length = [lx, 4 * radius, 4 * radius]
        cylinder_p = [-100, 2 * radius, 2 * radius]

    elif periodic == "y":
        lx = (boundaries["x"][0] + radius) - (boundaries["x"][1] - radius)
        ly = cell["y"][2]
        lz = (boundaries["z"][0] + radius) - (boundaries["z"][1] - radius)
        shift["x"] = boundaries["x"][1] - radius
        shift["z"] = boundaries["z"][1] - radius
        dx = 0
        dy = 1
        dz = 0
        print(coord, shift)
        # cx = coord["x"] - shift["x"]
        # cy = -ly
        # cz = coord["z"] - shift["z"]
        lc = ly
        length = [4 * radius, ly, 4 * radius]
        cylinder_p = [2 * radius, -100, 2 * radius]

    else:
        lx = (boundaries["x"][0] + radius) - (boundaries["x"][1] - radius)
        ly = (boundaries["y"][0] + radius) - (boundaries["y"][1] - radius)
        lz = cell["z"][2]
        shift["x"] = boundaries["x"][1] - radius
        shift["y"] = boundaries["y"][1] - radius
        dx = 0
        dy = 0
        dz = 1
        # cx = coord["x"] - shift["x"]
        # cy = coord["y"] - shift["y"]
        # cz = -lz
        lc = lz
        length = [4 * radius, 4 * radius, lz]
        cylinder_p = [2 * radius, 2 * radius, -100]

    cylinder_d = [dx, dy, dz]

    print("Polymod cell: ", length)  # , 'Shift cell: ', shift)
    return shift, length, cylinder_p, cylinder_d, lc


def run_data4lammps(charge, ffname, cell):

    directory = "./data4lammps/"
    path_r = os.path.join(directory, "main.py")
    path_2 = os.path.join(directory, "doAtomTyping.py")
    new_typing_command = "python3.6 {0} {1}".format(path_2, ffname)
    data4lammps_command = "python {0} {1} {2} {3} {4} {5} {6} {7}".format(
        path_r, 0, cell[0], 0, cell[1], 0, cell[2], charge
    )
    return new_typing_command, data4lammps_command


def read_dump(ifile, multiframe=0):
    Dir = {}
    src = open(ifile)
    current = ""
    timestep = 0
    for line in src.readlines():
        ln = line.split()
        if ln[0] == "ITEM:":
            if ln[1] == "TIMESTEP" or ln[1] == "BOX" or ln[1] == "ATOMS":
                current = ln[1]
                if ln[1] == "BOX":
                    Dir[timestep][current] = []
                if ln[1] == "ATOMS":
                    Dir[timestep][current] = {}
                    Dir[timestep][current]["id"] = 0
                    dump = {}
                    for i, j in enumerate(ln[2:]):
                        dump[i] = j
            if ln[1:] == "NUMBER OF ATOMS".split():
                current = "NUMBER OF ATOMS"
            continue
        if current == "TIMESTEP":
            timestep = eval(ln[0])
            if timestep not in Dir:
                Dir[timestep] = {}
        if current == "NUMBER OF ATOMS":
            Dir[timestep][current] = eval(ln[0])
        if current == "BOX":
            Dir[timestep][current] += [eval(i) for i in ln]
        if current == "ATOMS":
            cid = Dir[timestep][current]["id"]
            Dir[timestep][current][cid] = {}
            for i, j in enumerate(
                [(eval(k) if k[0].isdigit or k[0] == "-" else k) for k in ln]
            ):
                Dir[timestep][current][cid][dump[i]] = j
            Dir[timestep][current]["id"] += 1
    return Dir


# def dump2data(datafile, dumpfile, ofile):
#    Dir = read_data(datafile)
#    if "Velocities" not in Dir:
#        Dir["Velocities"] = {}
#    src = open(dumpfile)
#    box = "x"
#    for line in src.readlines():
#        ln = line.split()
#        if len(ln) == 2 and ln[0] != "ITEM:":
#            Dir[box + "lo"] = eval(ln[0])
#            Dir[box + "hi"] = eval(ln[1])
#            if box == "x":
#                box = "y"
#                continue
#            if box == "y":
#                box = "z"
#                continue
#            continue
#        if len(ln) == 8:
#            Dir["Atoms"][eval(ln[0])][3] = eval(ln[2])
#            Dir["Atoms"][eval(ln[0])][4] = eval(ln[3])
#            Dir["Atoms"][eval(ln[0])][5] = eval(ln[4])
#            Dir["Atoms"][eval(ln[0])][6] = 0
#            Dir["Atoms"][eval(ln[0])][7] = 0
#            Dir["Atoms"][eval(ln[0])][8] = 0

#      if eval(ln[0]) not in Dir['Velocities']:
#        Dir['Velocities'][eval(ln[0])]=[]
#        Dir['Velocities'][eval(ln[0])].append(eval(ln[5]))
#        Dir['Velocities'][eval(ln[0])].append(eval(ln[6]))
#        Dir['Velocities'][eval(ln[0])].append(eval(ln[7]))
#      Dir['Velocities'][eval(ln[0])][0]=eval(ln[5])
#      Dir['Velocities'][eval(ln[0])][1]=eval(ln[6])
#      Dir['Velocities'][eval(ln[0])][2]=eval(ln[7])
#    if len(ln)==11:
#      Dir['Atoms'][eval(ln[0])][3]=eval(ln[2])
#      Dir['Atoms'][eval(ln[0])][4]=eval(ln[3])
#      Dir['Atoms'][eval(ln[0])][5]=eval(ln[4])
#      Dir['Atoms'][eval(ln[0])][6]=eval(ln[8])
#      Dir['Atoms'][eval(ln[0])][7]=eval(ln[9])
#      Dir['Atoms'][eval(ln[0])][8]=eval(ln[10])
#      if eval(ln[0]) not in Dir['Velocities']:
#        Dir['Velocities'][eval(ln[0])]=[]
#        Dir['Velocities'][eval(ln[0])].append(eval(ln[5]))
#        Dir['Velocities'][eval(ln[0])].append(eval(ln[6]))
#        Dir['Velocities'][eval(ln[0])].append(eval(ln[7]))
#      Dir['Velocities'][eval(ln[0])][0]=eval(ln[5])
#      Dir['Velocities'][eval(ln[0])][1]=eval(ln[6])
#      Dir['Velocities'][eval(ln[0])][2]=eval(ln[7])
#    src.close()
#    write_data(Dir, ofile)


def write_data(Dir, ofile, v=1, a=0):
    des = open(ofile, "w")
    des.write("LAMMPS data file via Tongtong\n")
    des.write("\n")
    if a:
        ilist = ["atom"]
        List = ["Masses", "Atoms", "Velocities"]
    else:
        ilist = ["atom", "bond", "angle", "dihedral", "improper"]
        List = [
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
    for i in ilist:
        if (i + "s") in Dir:
            des.write("%d %s\n" % (Dir[i + "s"], (i + "s")))
            des.write("%d %s\n" % (Dir[i + " types"], (i + " types")))
    des.write("\n")
    for i in ["x", "y", "z"]:
        des.write(
            "%f %f %s %s\n" % (Dir[i + "lo"], Dir[i + "hi"], (i + "lo"), (i + "hi"))
        )
    des.write("\n")
    if v == 0:
        List.remove("Velocities")
    for key in List:
        if key in Dir and len(Dir[key]) > 0:
            des.write(key + "\n")
            des.write("\n")
            for i in Dir[key]:
                des.write(str(i) + " " + " ".join(str(j) for j in Dir[key][i]) + "\n")
            des.write("\n")
    des.close()


def input_coat4(na, cpos, idata, odata):
    des = open("run_coat4.in", "w")
    des.write("# General parameters\n")
    des.write("units real\n")
    des.write("atom_style        full\n")
    des.write("boundary          p p p\n")
    des.write("special_bonds     lj/coul 0.0 0.0 1.0 dihedral yes\n")
    des.write("dielectric        1.0\n")
    des.write("pair_style        lj/cut/coul/long  12.0\n")
    des.write("bond_style        harmonic\n")
    des.write("angle_style       harmonic\n")
    des.write("dihedral_style    harmonic\n")
    des.write("improper_style    harmonic\n")
    des.write("kspace_style      pppm 1.0e-6\n")
    des.write("\n")
    des.write("read_data         %s\n" % idata)
    des.write("\n")
    des.write("thermo_style      custom  temp vol density pxx pyy pzz lx ly lz\n")
    des.write("thermo            100\n")
    des.write("thermo_modify     flush yes line multi\n")
    des.write("\n")
    des.write("dump              1 all atom 500 coat2.dump\n")
    des.write("\n")
    des.write("label       loopa\n")
    des.write("variable    i loop %d\n" % na)
    des.write("  variable  k equal $i\n")
    des.write("  variable  h equal x[$k]\n")
    des.write("  variable  p equal y[$k]\n")
    des.write('  print     "x= $h"\n')
    des.write('  print     "y= $p"\n')
    des.write("  variable  dx equal $h-%f\n" % cpos[0])
    des.write("  variable  dy equal $p-%f\n" % cpos[1])
    des.write("  variable  dr equal sqrt(${dx}*${dx}+${dy}*${dy})\n")
    des.write("  variable  nvx equal -${dx}/${dr}*0.01\n")
    des.write("  variable  nvy equal -${dy}/${dr}*0.01\n")
    des.write("  set atom $k vx ${nvx} vy ${nvy} vz 0.0\n")
    des.write("next        i\n")
    des.write("jump        SELF loopa\n")
    des.write("\n")
    des.write("write_data        %s\n" % odata)
    des.write("\n")
    des.close()

    # with open('coat4.in', 'r') as file:
    # filedata = file.read()

    # Replace the target string
    # filedata = filedata.replace('na', str(na))
    # filedata = filedata.replace('ccx', str(cpos[0]))
    # filedata = filedata.replace('ccy', str(cpos[1]))
    # filedata = filedata.replace('idata', idata)
    # filedata = filedata.replace('odata', odata)

    # Write the file out again
    # with open('run_coat4.in', 'w') as file:
    # file.write(filedata)


def input_coat5(cpos, radius, oradius, idata, odata, tstep, X6paircoeffs, mctemp):
    des = open("run_coat5.in", "w")
    des.write("# General parameters\n")
    des.write("units  real\n")
    des.write("atom_style        full\n")
    des.write("boundary          p p p\n")
    des.write("special_bonds     lj/coul 0.0 0.0 1.0 dihedral yes\n")
    des.write("dielectric        1.0\n")
    des.write("#pair_style        lj/cut/coul/long  12.0\n")
    des.write("pair_style       buck/coul/cut 10.0\n")
    des.write("bond_style        harmonic\n")
    des.write("angle_style       harmonic\n")
    des.write("dihedral_style    harmonic\n")
    des.write("improper_style    harmonic\n")
    des.write("\n")
    des.write("read_data         %s\n" % idata)
    des.write("\n")
    src = open(X6paircoeffs)
    for line in src.readlines():
        des.write(line)
    src.close()
    des.write("\n")
    des.write("thermo_style      custom  temp vol density pxx pyy pzz lx ly lz\n")
    des.write("thermo            100\n")
    des.write("thermo_modify     flush yes line one\n")
    des.write("\n")
    des.write(
        "region inner cylinder z %f %f %f INF INF units box side out\n"
        % (cpos[0], cpos[1], radius - 2)
    )
    des.write(
        "region outter cylinder z %f %f %f INF INF units box side in\n"
        % (cpos[0], cpos[1], oradius + 3.5)
    )
    des.write("\n")
    des.write("comm_style tiled\n")
    des.write("fix LB all balance 1000 1.1 rcb\n")
    des.write("\n")
    des.write("fix               1 all wall/region inner lj126 5.0 3.5 12.0\n")
    des.write("fix               8 all wall/region outter lj126 5.0 3.5 3.5\n")
    des.write("\n")
    des.write("fix               2 all nvt temp %f %f 100 \n" % (mctemp, mctemp))
    des.write("\n")
    des.write("reset_timestep %d\n" % tstep)
    des.write("\n")
    des.write("dump              1 all atom 500 coat5.*.dump\n")
    des.write("\n")
    des.write("run 5000\n")
    des.write("write_data        %s\n" % odata)
    des.write("\n")
    des.close()

    # with open('coat5.in', 'r') as file:
    # filedata = file.read()

    # Replace the target string
    # filedata = filedata.replace('ccx', str(cpos[0]))
    # filedata = filedata.replace('ccy', str(cpos[1]))
    # filedata = filedata.replace('radius', str(radius-2))
    # filedata = filedata.replace('oadius', str(oradius+3.5))
    # filedata = filedata.replace('idata', idata)
    # filedata = filedata.replace('odata', odata)
    # filedata = filedata.replace('tstep', str(tstep))

    # Write the file out again
    # with open('run_coat5.in', 'w') as file:
    # file.write(filedata)


def input_coat6(na, idata, odata):
    des = open("run_coat6.in", "w")
    des.write("# General parameters\n")
    des.write("units  real\n")
    des.write("atom_style        full\n")
    des.write("boundary          p p p\n")
    des.write("#special_bonds     lj/coul 0.0 0.0 1.0 dihedral yes\n")
    des.write("#dielectric        1.0\n")
    des.write("#pair_style        lj/cut/coul/long  12.0\n")
    des.write("#bond_style        harmonic\n")
    des.write("#angle_style       harmonic\n")
    des.write("#improper_style    harmonic\n")
    des.write("#kspace_style      pppm 1.0e-6\n")
    des.write("\n")
    des.write("read_data %s\n" % idata)
    des.write("\n")
    des.write("thermo_style      custom  temp vol density pxx pyy pzz lx ly lz\n")
    des.write("thermo            100\n")
    des.write("thermo_modify     flush yes line multi\n")
    des.write("\n")
    des.write("dump              1 all atom 500 coat2.dump\n")
    des.write("\n")
    des.write("label       loopa\n")
    des.write("variable    i loop %d\n" % na)
    des.write("  variable  k equal $i\n")
    des.write("  variable  h equal x[$k]\n")
    des.write("  variable  p equal y[$k]\n")
    des.write('  print     "x= $h"\n')
    des.write('  print     "y= $p"\n')
    des.write("  variable  dx equal $h \n")
    des.write("  variable  dy equal $p \n")
    des.write("  variable  dr equal sqrt(${dx}*${dx}+${dy}*${dy})\n")
    des.write("  variable  nvx equal -${dx}/${dr}*0.001\n")
    des.write("  variable  nvy equal -${dy}/${dr}*0.001\n")
    des.write("  set atom $k vx ${nvx} vy ${nvy} vz 0.0\n")
    des.write("next        i\n")
    des.write("jump        SELF loopa\n")
    des.write("\n")
    des.write("write_data        %s\n" % odata)
    des.write("\n")
    des.close()

    # with open('coat6.in', 'r') as file:
    # filedata = file.read()

    # Replace the target string
    # filedata = filedata.replace('na', str(na))
    # filedata = filedata.replace('idata', idata)
    # filedata = filedata.replace('odata', odata)

    # Write the file out again
    # with open('run_coat6.in', 'w') as file:
    # file.write(filedata)


def input_coat7(oradius, idata, odata, tstep, polymer_type, HE_type, mctemp):
    des = open("run_coat7.in", "w")
    des.write("newton          on\n")
    des.write("boundary        p p p\n")
    des.write("units           real\n")
    des.write("box tilt large\n")
    des.write("\n")
    des.write("include      ../potential_head.mod\n")
    des.write("read_data  %s\n" % idata)
    des.write("\n")
    des.write("group    polymer type 1:%d\n" % polymer_type)
    des.write("group    HE type %d:%d\n" % (polymer_type + 1, HE_type + polymer_type))
    des.write("\n")
    des.write("include      potential.mod\n")
    des.write("\n")
    des.write("  compute                              stress  all  stress/atom NULL\n")
    des.write("  compute                              PEbond all  pe/atom bond\n")
    des.write("  compute                              PEangle all  pe/atom angle\n")
    des.write(
        "  compute                              PEdihed               all  pe/atom dihedral\n"
    )
    des.write("  compute                              PEimp   all  pe/atom improper\n")
    des.write("  compute                              PEinter all  pe/atom pair\n")
    des.write("\n")
    des.write(
        "thermo_style     custom step time etotal pe ke temp press pxx pyy pzz pxy pxz pyz density evdwl ecoul epair \
            ebond eangle edihed eimp lx ly evdwl\n"
    )
    des.write("thermo                 5\n")
    des.write("\n")
    des.write("#-------------------------------------------------------------\n")
    des.write("# SIMULATION SETUP\n")
    des.write("#-------------------------------------------------------------\n")
    des.write("#\n")
    des.write("\n")
    des.write("comm_style tiled\n")
    des.write("fix 2 all balance 1000 1.1 rcb\n")
    des.write("\n")
    des.write("run_style  verlet\n")
    des.write("\n")
    des.write(
        "region outter cylinder z 0.0 0.0 %f INF INF units box side in\n"
        % (oradius + 3.5)
    )
    des.write("\n")
    des.write(
        "fix              6 HE rigid group 1 HE force * off off off torque * off off off\n"
    )
    des.write("fix    7 polymer nvt temp %f %f 50\n" % (mctemp, mctemp))
    des.write("fix               8 all wall/region outter lj126 5.0 3.5 3.5\n")
    des.write("dump              1 all atom 500 coat7.*.dump\n")
    des.write("\n")
    des.write("reset_timestep %d\n" % tstep)
    des.write("\n")
    des.write("run             5000\n")
    des.write("unfix           6\n")
    des.write("unfix           7\n")
    des.write("\n")
    des.write("write_data %s\n" % odata)
    des.write("\n")
    des.close()

    # with open('coat7.in', 'r') as file:
    # filedata = file.read()

    # Replace the target string
    # filedata = filedata.replace('oadius', str(oradius+3.5))
    # filedata = filedata.replace('idata', idata)
    # filedata = filedata.replace('odata', odata)
    # filedata = filedata.replace('tstep', str(tstep))
    # filedata = filedata.replace('polymer_type', ('1:%d' %polymer_type))
    # filedata = filedata.replace('HE_type', ('%d:%d' %(polymer_type+1,HE_type+polymer_type)))

    # Write the file out again
    # with open('run_coat7.in', 'w') as file:
    # file.write(filedata)


def Get_Mass_Radius(Dir, center, dimension=3, plane="xy"):
    D = 0
    if dimension == 3:
        for a_id in Dir["Atoms"]:
            d = 0
            for i in range(3):
                d += (Dir["Atoms"][a_id][3 + i] - center[i]) ** 2
            if D < d:
                D = d
    else:
        i_range = []
        index = {"x": 0, "y": 1, "z": 2}
        for xyz in ["x", "y", "z"]:
            if xyz in plane:
                i_range.append(index[xyz])
        for a_id in Dir["Atoms"]:
            d = 0
            for i in i_range:
                d += (Dir["Atoms"][a_id][3 + i] - center[i]) ** 2
            if D < d:
                D = d
    return np.sqrt(D), center


def Get_Center_of_Mass(Dir_ori):
    M = 0
    Mx = 0
    My = 0
    Mz = 0
    for a_id in Dir_ori["Atoms"]:
        m = Dir_ori["Masses"][Dir_ori["Atoms"][a_id][1]][0]
        M += m
        Mx += m * Dir_ori["Atoms"][a_id][3]
        My += m * Dir_ori["Atoms"][a_id][4]
        Mz += m * Dir_ori["Atoms"][a_id][5]
    return [Mx / M, My / M, Mz / M]


def data_Translation(Dir, vector, box=0):
    for a_id in Dir["Atoms"]:
        Dir["Atoms"][a_id][3] += vector[0]
        Dir["Atoms"][a_id][4] += vector[1]
        Dir["Atoms"][a_id][5] += vector[2]
    if box:
        for i, xyz in enumerate(["x", "y", "z"]):
            for lh in ["lo", "hi"]:
                Dir[xyz + lh] += vector[i]
    return Dir


def add_data(Dir_1, Dir_2, add="append"):
    Dir_data = {}
    for i in ["atom", "bond", "angle", "dihedral", "improper"]:
        if (i + "s") in Dir_1:
            Dir_data[i + "s"] = Dir_1[i + "s"]
            Dir_data[i + " types"] = Dir_1[i + " types"]
        if (i + "s") in Dir_2:
            if (i + "s") in Dir_data:
                Dir_data[i + "s"] += Dir_2[i + "s"]
                Dir_data[i + " types"] += Dir_2[i + " types"]
            else:
                Dir_data[i + "s"] = Dir_2[i + "s"]
                Dir_data[i + " types"] = Dir_2[i + " types"]
        if ("extra " + i + " per atom") in Dir_1:
            Dir_data["extra " + i + " per atom"] = Dir_1["extra " + i + " per atom"]
        if ("extra " + i + " per atom") in Dir_2:
            if ("extra " + i + " per atom") in Dir_data:
                Dir_data["extra " + i + " per atom"] = max(
                    Dir_1["extra " + i + " per atom"], Dir_2["extra " + i + " per atom"]
                )
            else:
                Dir_data["extra " + i + " per atom"] = Dir_2["extra " + i + " per atom"]
    for i in ["x", "y", "z"]:
        if (i + "lo") in Dir_1:
            Dir_data[i + "lo"] = Dir_1[i + "lo"]
            Dir_data[i + "hi"] = Dir_1[i + "hi"]
        if (i + "lo") in Dir_2:
            if (i + "lo") not in Dir_data:
                Dir_data[i + "lo"] = Dir_2[i + "lo"]
                Dir_data[i + "hi"] = Dir_2[i + "hi"]
    List = [
        "Masses",
        "Atoms",
        "Velocities",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
    ]
    for key in List:
        if key in Dir_1 and len(Dir_1[key]) > 0:
            Dir_data[key] = Dir_1[key]
        if key in Dir_2 and len(Dir_2[key]) > 0:
            if key in Dir_data:
                if key in [
                    "Masses",
                    "Atoms",
                    "Velocities",
                    "Bonds",
                    "Angles",
                    "Dihedrals",
                    "Impropers",
                ]:
                    if key == "Atoms":
                        a_ids = {}
                        inimol = max([Dir_data[key][i][0] for i in Dir_data[key]])
                    ini = max([i for i in Dir_data[key]])
                    for a_id in Dir_2[key]:
                        if key == "Atoms":
                            a_ids[a_id] = ini + a_id
                        Dir_data[key][ini + a_id] = Dir_2[key][a_id]
                        if key == "Atoms":
                            Dir_data[key][ini + a_id][0] = Dir_2[key][a_id][0] + inimol
                            Dir_data[key][ini + a_id][1] = (
                                Dir_2[key][a_id][1] + Dir_1["atom types"]
                            )
                        if key == "Bonds":
                            Dir_data[key][ini + a_id][0] = (
                                Dir_2[key][a_id][0] + Dir_1["bond types"]
                            )
                        if key == "Angles":
                            Dir_data[key][ini + a_id][0] = (
                                Dir_2[key][a_id][0] + Dir_1["angle types"]
                            )
                        if key == "Dihedrals":
                            Dir_data[key][ini + a_id][0] = (
                                Dir_2[key][a_id][0] + Dir_1["dihedral types"]
                            )
                        if key == "Impropers":
                            Dir_data[key][ini + a_id][0] = (
                                Dir_2[key][a_id][0] + Dir_1["improper types"]
                            )
                        if key in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
                            for i in range(1, len(Dir_2[key][a_id])):
                                Dir_data[key][ini + a_id][i] = a_ids[
                                    Dir_data[key][ini + a_id][i]
                                ]
            else:
                Dir_data[key] = Dir_2[key]
    return Dir_data


# def grabprocessors(ifile):
#    src = open(ifile)
#    for line in src.readlines():
#        ln = line.split()
#        if ln and ln[0] == "read_data":
#            datafile = ln[1]
#            break
#    src.close()
#    Dir = read_data(datafile)
#    atoms = Dir["atoms"]
#    if "bonds" in Dir:
#        return atoms / 2000 + 1
#    else:
#        return atoms / 1000 + 1


def correctppn(ppn):
    nodes = ppn / 20 + 1 if (ppn / 20 and ppn % 20) else ppn / 20 if ppn / 20 else 1
    if not ppn / 20:
        ppn = (
            20
            if ppn > 10
            else 10
            if ppn > 8
            else 8
            if ppn > 4
            else 4
            if ppn > 2
            else 2
            if ppn > 1
            else 1
        )
    return nodes, 20 if ppn / 20 else ppn, nodes * 20 if ppn / 20 else ppn


def shiftpotential(ifile1, ifile2, ofile, Dir_polymer):
    src = open(ifile1)
    des = open(ofile, "w")
    for line in src.readlines():
        ln = line.split("#")[0].split()
        if ln:
            if ln[0] == "pair_coeff":
                des.write(" ".join(ln[0:3]) + " buck " + " ".join(ln[3:]) + "\n")
    src.close()
    src = open(ifile2)
    for line in src.readlines():
        ln = line.split("#")[0].split()
        if ln:
            if ln[0] == "pair_coeff":
                if ln[1] != "*":
                    ln[1] = str(eval(ln[1]) + Dir_polymer["atom types"])
                if ln[2] != "*":
                    ln[2] = str(eval(ln[2]) + Dir_polymer["atom types"])
            if ln[0] == "bond_coeff":
                ln[1] = str(eval(ln[1]) + Dir_polymer["bond types"])
            if ln[0] == "angle_coeff":
                ln[1] = str(eval(ln[1]) + Dir_polymer["angle types"])
            if ln[0] == "dihedral_coeff":
                ln[1] = str(eval(ln[1]) + Dir_polymer["dihedral types"])
            if ln[0] == "improper_coeff":
                ln[1] = str(eval(ln[1]) + Dir_polymer["improper types"])
        des.write(" ".join(ln) + "\n")
    src.close()
    for key in ["Bond", "Angle", "Dihedral", "Improper"]:
        if (key + " Coeffs") in Dir_polymer:
            for i in Dir_polymer[key + " Coeffs"]:
                if key == "Bond":
                    des.write(
                        "bond_coeff %i harmonic %s\n"
                        % (
                            i,
                            " ".join([str(k) for k in Dir_polymer[key + " Coeffs"][i]]),
                        )
                    )
                elif key == "Angle":
                    des.write(
                        "angle_coeff %i %s\n"
                        % (
                            i,
                            " ".join([str(k) for k in Dir_polymer[key + " Coeffs"][i]]),
                        )
                    )
                elif key == "Dihedral":
                    des.write(
                        "dihedral_coeff %i %s\n"
                        % (
                            i,
                            " ".join([str(k) for k in Dir_polymer[key + " Coeffs"][i]]),
                        )
                    )
                elif key == "Improper":
                    des.write(
                        "improper_coeff %i %s\n"
                        % (
                            i,
                            " ".join([str(k) for k in Dir_polymer[key + " Coeffs"][i]]),
                        )
                    )
    des.close()


def addatomtype(ifile, ofile, elementlist):
    src = open(ifile)
    des = open(ofile, "w")
    a = 1
    for line in src.readlines():
        ln = line.split()
        if ln:
            des.write(" ".join(ln) + "\n")
            a += 1
    src.close()
    for e in elementlist:
        if e == "C":
            des.write("%d C_3\n" % a)
        else:
            des.write("%d %s\n" % (a, e))
        a += 1
    des.close()


# def writeghostdata(ifile):
#    Dir = read_data(ifile)
#    write_data(Dir, "ghost_grain.data", v=0, a=1)


def writepolymodfile(polymer_type_custom, ofile1):
    mw = {}
    des = open(ofile1, "w")
    des.write("#\n")
    des.write("# PolymerBuilder input file: order of inputs is critical\n")
    des.write("#\n")
    des.write("\n")
    des.write("# Temperature (K)\n")
    des.write("mctemp  # Temperature (K)\n")
    des.write("\n")
    des.write("# Elements: atomic type indices for all monomers\n")
    des.write(
        "%d                    # Number of atom types\n"
        % polymer_type_custom["atom_type"]
    )
    atom_type = {}
    a_id = 0
    for i in range(polymer_type_custom["nm"]):
        for j in range(len(polymer_type_custom["monomer"][i + 1]["info"])):
            if polymer_type_custom["monomer"][i + 1]["info"][j][0] not in atom_type:
                a_id += 1
                atom_type[a_id] = polymer_type_custom["monomer"][i + 1]["info"][j][0]
                atom_type[polymer_type_custom["monomer"][i + 1]["info"][j][0]] = a_id
    ofile = open("str2lammps/types/custom_atom_type.dat", "w")
    for i in range(polymer_type_custom["atom_type"]):
        des.write(
            "%s                    # Type %d element\n"
            % (atom_type[i + 1].split("_")[0], (i + 1))
        )
        ofile.write("%d %s\n" % ((i + 1), atom_type[i + 1]))
    ofile.close()
    des.write("\n")
    des.write("# Number of monomers\n")
    des.write("%d  # Number of monomer types\n" % polymer_type_custom["nm"])
    des.write("\n")
    for i in range(polymer_type_custom["nm"]):
        des.write("# Monomer: z-matrix\n")
        des.write(
            "%d                    # Number of atoms in z-matrix (lengths in A, angles in degrees)\n"
            % len(polymer_type_custom["monomer"][i + 1]["info"])
        )
        for j in range(len(polymer_type_custom["monomer"][i + 1]["info"])):
            des.write(
                str(atom_type[polymer_type_custom["monomer"][i + 1]["info"][j][0]])
                + " "
                + " ".join(polymer_type_custom["monomer"][i + 1]["info"][j][1:])
                + "\n"
            )
        des.write(
            "%d                    # Number of backbone atoms, z-matrix top\n"
            % polymer_type_custom["monomer"][i + 1]["torsion"]["len"]
        )
        des.write("\n")
        des.write("# Monomer: backbone torsion angle probabilities\n")
        des.write(
            "# Starting with backbone atom 3, specify whether the torsion should change, \n"
        )
        des.write("#   and, if so, how.\n")
        des.write("# Values for specification:\n")
        des.write("#   0: no change\n")
        des.write("#   1: uniform probability for all angles\n")
        des.write("#   2: energy levels associated with specific angles\n")
        des.write("#   3: probability associated with specific angles\n")
        des.write("\n")
        if "all" in polymer_type_custom["monomer"][i + 1]["torsion"]:
            for j in range(polymer_type_custom["monomer"][i + 1]["torsion"]["len"] - 2):
                tor_type = int(
                    polymer_type_custom["monomer"][i + 1]["torsion"]["all"][0]
                )
                des.write(
                    "%d                    # Torsion %d specification\n"
                    % (tor_type, (j + 3))
                )
                if tor_type == 2 or tor_type == 3:
                    src = open(
                        "../"
                        + polymer_type_custom["monomer"][i + 1]["torsion"]["all"][1]
                    )
                    for line in src.readlines():
                        des.write(line)
                    src.close()
                    des.write("\n")
        else:
            for j in range(polymer_type_custom["monomer"][i + 1]["torsion"]["len"] - 2):
                if str(j + 3) in polymer_type_custom["monomer"][i + 1]["torsion"]:
                    tor_type = int(
                        polymer_type_custom["monomer"][i + 1]["torsion"][str(j + 3)][0]
                    )
                    des.write(
                        "%d                    # Torsion %d specification\n"
                        % (tor_type, (j + 3))
                    )
                    if tor_type == 2 or tor_type == 3:
                        src = open(
                            "../"
                            + polymer_type_custom["monomer"][i + 1]["torsion"][
                                str(j + 3)
                            ][1]
                        )
                        for line in src.readlines():
                            des.write(line)
                        src.close()
                        des.write("\n")
                else:
                    print(
                        "Backbone torsion angle %d probabilities type for monomer %d not specified, use default with \
                            no change"
                        % ((j + 3), (i + 1))
                    )
                    tor_type = 0
                    des.write(
                        "%d                    # Torsion %d specification: no change\n"
                        % (tor_type, (j + 3))
                    )
        des.write("\n")
    des.write(
        "3.0                    # Torsion delta: change in torsions is +/- this value\n"
    )
    des.write(
        "10                   # Number of torsion delta steps to minimize torsion energy\n"
    )
    des.write("\n")
    des.write("# Backbone bond length between all monomers (A)\n")
    des.write("1.53\n")
    des.write("\n")
    des.write("# Monomer arrangements\n")
    des.write(
        "%d                    # Number of monomer arrangements\n"
        % polymer_type_custom["nc"]
    )
    for i in range(polymer_type_custom["nc"]):
        mw[i] = 0
        if polymer_type_custom["chain"][i + 1]["arrangement"]["type"]:
            des.write(
                "1                    # Arrangement: 0 = pattern, 1 = probability\n"
            )
            des.write(
                "%s              # Probability of monomer(s)\n"
                % (
                    " ".join(
                        polymer_type_custom["chain"][i + 1]["arrangement"]["sequence"]
                    )
                )
            )
            length = len(polymer_type_custom["chain"][i + 1]["arrangement"]["sequence"])
            for j in range(length):
                mw[i] += polymer_type_custom["monomer"][j + 1]["mass"] * float(
                    polymer_type_custom["chain"][i + 1]["arrangement"]["sequence"][j]
                )
        else:
            des.write(
                "0                    # Arrangement: 0 = pattern, 1 = probability\n"
            )
            des.write(
                "%d                    # Number of monomers in first pattern\n"
                % polymer_type_custom["chain"][i + 1]["arrangement"]["len"]
            )
            des.write(
                "%s              # Repeat...\n"
                % (
                    " ".join(
                        polymer_type_custom["chain"][i + 1]["arrangement"]["sequence"]
                    )
                )
            )
            length = polymer_type_custom["chain"][i + 1]["arrangement"]["len"]
            for j in polymer_type_custom["chain"][i + 1]["arrangement"]["sequence"]:
                mw[i] += polymer_type_custom["monomer"][int(j)]["mass"] / length
    r_mw = 0
    for i in range(polymer_type_custom["nc"]):
        des.write("%f " % polymer_type_custom["chain"][i + 1]["probability"])
        r_mw += mw[i] * polymer_type_custom["chain"][i + 1]["probability"]
    des.write("             # Probabilities of each monomer arrangement\n")
    des.write("\n")
    des.write("# System\n")
    des.write("nc                   # Number of chains to build\n")
    des.write("nm                   # Number of monomers per chain\n")
    des.write(
        "lcx lcy lcz             # Dimensions of cell (A); not used if density > 0.0\n"
    )
    des.write("0.0                 # Density in g/cm^3\n")
    des.write("\n")
    des.write("# Excluded cylinders -- define volumes in which no polymer exists\n")
    des.write("1                    # Number of excluded cylinders\n")
    des.write(
        "ccx ccy ccz             # Cylinder 1: start position x y z -- cylinder end center \n"
    )
    des.write(
        "dx dy dz              # Cylinder 1: direction from start position -- axis\n"
    )
    des.write("crad                  # Cylinder 1: radius, extension from axis\n")
    des.write("clc                 # Cylinder 1: length\n")
    des.write("\n")
    des.write("1 # Included cylinders -- define volumes in which polymers exist\n")
    des.write(
        "ccx ccy ccz   # Cylinder 1 : start position x y z -- cylinder end center\n"
    )
    des.write("dx dy dz   # Cylinder 1 : direction from start position -- axis\n")
    des.write("irad              # Cylinder 1 : radius, extension from axis\n")
    des.write("clc              # Cylinder 1 : length\n")
    des.write("\n")
    des.write("# Excluded slabs -- define volumes in which no polymer exists\n")
    des.write("0                    # Number of excluded slabs\n")
    des.write("\n")
    des.write("# Configurations and interactions\n")
    des.write("40                   # Max number of configurations to test\n")
    des.write(
        "1.0                  # Keep chains which are this fraction of desired length\n"
    )
    des.write("1                    # 1: Self-avoiding\n")
    des.write("1                  # Self-avoiding cutoff (A)\n")
    des.write("1                    # 1: Long range interactions\n")
    des.write("5                  # Interaction cutoff (A)\n")
    des.write("5.0                  # Bin size (A)\n")
    des.write("4                    # Bond cutoff\n")
    des.write("\n")
    des.write("# Output\n")
    des.write(
        "1                    # 1: Write (unwrapped) PDB file of constructed chains\n"
    )
    des.write(
        "1                    # 1: Write wrapped PDB file of constructed chains\n"
    )
    des.write("0                    # 1: Write PDB files of monomer rotation\n")
    des.write(
        "0                    # 1: Write output file of chain length vs. # monomers\n"
    )
    des.write("0                    # 1: Write output file of chain length histogram\n")
    des.write(
        "0                    # 1: Write torsion angle probabilities, selection histogram\n"
    )
    des.write("0                    # 1: Write z-matrices of all chains\n")
    des.write("0                    # 1: Write final system energy\n")
    des.write("\n")
    des.write("#\n")
    des.write("# Status and messages\n")
    des.write("#\n")
    des.write("1                    # 1: Write messages to stdout;  0: Write to file\n")
    des.write(
        "# If previous line is 1, nothing more is needed; if it is 0, file name follows\n"
    )
    des.write("#/path/to/log.build\n")
    des.write("1                    # 1: Write status to stdout;    0: Write to file\n")
    des.write(
        "# XXX Only one of status.build or bar/Rappture flag should be uncommented!!\n"
    )
    des.write("# If writing status to file, file name follows\n")
    des.write("#/path/to/status.build\n")
    des.write(
        "# ELSE if writing status to stdout, 1: terminal bar, 0: Rappture status lines\n"
    )
    des.write("1\n")
    des.write("\n")
    des.write("# RNG seed; use time if 0\n")
    des.write("0                    \n")
    des.write("\n")
    des.write(
        "# Scale factor used to compare atom distances in monomer to equilibrium bond\n"
    )
    des.write(
        "# distances when searching for bonds not specified in monomer z-matrix\n"
    )
    des.write("1.1\n")
    des.write("\n")
    des.close()
    return r_mw


def write_minimize(ofile, X6file, afile):
    des = open(ofile, "w")
    des.write("# General parameters\n")
    des.write("units  real\n")
    des.write("atom_style        full\n")
    des.write("boundary          p p p\n")
    des.write("special_bonds     lj/coul 0.0 0.0 1.0 dihedral yes\n")
    des.write("dielectric        1.0\n")
    des.write("pair_style        lj/cut  12.0\n")
    des.write("bond_style        harmonic\n")
    des.write("angle_style       harmonic\n")
    des.write("dihedral_style    harmonic\n")
    des.write("improper_style    harmonic\n")
    des.write("read_data step0.data #        polymer_relax.data\n")
    des.write("neighbor          0.3 bin\n")
    des.write(
        "thermo_style      custom step etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press pxx pyy pzz \
        pxy pxz pyz lx ly lz vol density\n"
    )
    des.write("thermo            10\n")
    des.write("thermo_modify     flush yes\n")
    des.write("fix LB all balance 1000 1.05 shift xy 10 1.05\n")
    des.write("# Minimization parameters\n")
    des.write("min_style         cg  # hftn\n")
    des.write("min_modify        dmax 0.02\n")
    des.write("min_modify        line quadratic  # backtrack\n")
    des.write("neigh_modify      every 1 delay 0\n")
    LJ_params = {
        "H": [3.195, 0.0152],
        "C": [3.8983, 0.0951],
        "N": [3.6621, 0.0774],
        "O": [3.4046, 0.0957],
        "F": [3.4720, 0.0725],
        "S": [4.0300, 0.3440],
        "Cl": [3.9503, 0.2833],
        "Si": [4.27, 0.31],
    }
    d0_to_epsilon = 1.0
    r0_to_sigma = 1.0 / (2.0 ** (1.0 / 6.0))
    lj_D0 = {}
    lj_R0 = {}
    src = open(afile)
    el_types = {}
    for line in src.readlines():
        ln = line.split()
        if ln:
            el_types[int(ln[0])] = ln[1].split("_")[0]
    ntypes = len(el_types)
    for i in range(ntypes):
        lj_R0[i] = {}
        lj_D0[i] = {}
        param_i = LJ_params[el_types[i + 1]]
        for j in range(i, ntypes):
            param_j = LJ_params[el_types[j + 1]]
            lj_R0[i][j] = np.sqrt(param_i[0] * param_j[0])
            lj_D0[i][j] = np.sqrt(param_i[1] * param_j[1])
    lammps_min_steps = 5000
    lammps_min_levels = 3
    lammps_min_init = 0.5
    dp = (1.0 - lammps_min_init) / (lammps_min_levels - 1)
    for level in range(lammps_min_levels):
        des.write("# Minimization %d\n" % (i + 1))
        for i in range(ntypes):
            for j in range(i, ntypes):
                des.write("pair_coeff        %d  %d  " % (i + 1, j + 1))
                d0 = lj_D0[i][j] * (lammps_min_init + dp * level)
                r0 = lj_R0[i][j] * (lammps_min_init + dp * level)
                des.write("%f  %f\n" % (d0 * d0_to_epsilon, r0 * r0_to_sigma))
        des.write("minimize          1.0e-9 1.0e-9 %d 100000\n" % lammps_min_steps)
    des.write("# Dump minimized system\n")
    des.write("dump              1 all atom 1 min.dump\n")
    des.write("dump_modify       1 image yes scale no\n")
    des.write("run               0\n")
    des.write("undump            1\n")
    des.write("# MD parameters\n")
    des.write("neigh_modify every 1 delay 5\n")
    des.write("pair_style        buck/coul/long  12.0 12.0\n")
    des.write("kspace_style      pppm 1e-4\n")
    src = open(X6file)
    for line in src.readlines():
        des.write(line)
    des.write("minimize          1.0e-9 1.0e-9 %d 100000\n" % lammps_min_steps)
    des.write("write_data step1.data\n")
    des.write("\n")


# def main(args):

#    args = complete_args(args)

#    os.system("cp polymer_types/" + polymer_type + "/* .")

# Get number of chains, molecules
#    p_fra = 0.05
#    if "p_fra" in args:
#        if len(args["p_fra"]) == 1:
#            p_fra = eval(args["p_fra"][0])
#            if p_fra < 0 or p_fra > 1:
#                print("Please input a valid number between 0 to 1")
#                return False
#            else:
#                print(
#                    "Molecular weight fraction of polymer coated: "
#                    + str(p_fra * 100)
#                    + "%"
#                )
#        else:
#            print("Please specify molecular weight fraction of polymer coated")
#            return False
#    else:
#        print("Default 5% molecular weight fraction of polymer coated")

#    nm = 40
#    if "nm" in args:
#        if len(args["nm"]) == 1:
#            nm = eval(args["nm"][0])
#            if nm < 1:
#                print("Please input a valid number equal or larger than 1")
#                return False
#            else:
#                print(str(int(nm)) + " monomers per chain")
#        else:
#            print("Please specify number of monomers per chain")
#            return False
#    else:
#        print("Default 40 monomers per chain")

#    mctemp = 600
#    if "mctemp" in args:
#        if len(args["mctemp"]) == 1:
#            mctemp = eval(args["mctemp"][0])
#            if mctemp < 0:
#                print("Please input a valid number equal or larger than 0")
#                return False
#            else:
#                print("Monte Carlo temperature: " + str(mctemp) + " K")
#        else:
#            print("Please specify Monte Carlo temperature")
#            return False
#    else:
#        print("Default Monte Carlo temperature: 600 K")

#    np_flag = 0
#    rnp_flag = 0
#    if "parallel" in args:
#        np_flag = 1
#        print("LAMMPS will run in parallel mode")
#        if len(args["parallel"]) > 0:
#            if args["parallel"] == ["np"]:
#                if len(args["parallel"]) > 1:
#                    np = int(eval(args["parallel"][1]))
#                    if np <= 0:
#                        print("Please input a valid number larger than 0")
#                        return False
#                    else:
#                        print("%d of processors will be in use" % np)
#                else:
#                    print("Please input a number of processors you want to use")
#                    return False
#            else:
#                rnp_flag = 1
#                print(
#                    "Will calculating recommended number of processors after initial polymer configuration generated"
#                )
#        else:
#            rnp_flag = 1
#            print(
#                "Will calculating recommended number of processors after initial polymer configuration generated"
#            )

#    print("Running Lammps, get information about grain")
# Run Lammps, get information about grain
# Center of cylinder from upper and lower coordinates
# input_file(datafile)
# os.chdir('..')
#    writeinlammps(datafile, potential_headfile, potentialfile)
#    run_lammps()
# os.system('mv inr.lammps code')
# os.system('mv log.lammps code')
# os.system('mv tmp.out code')
# os.chdir('code')
#    cell_sizes, delta_cell = read_cell_sizes("../" + datafile)
#    periodic = min(delta_cell.items(), key=lambda x: x[1])[0]
#    mass, com, boundaries, coord, radius = get_boundaries(periodic)
#    dimension = 2
#    plane = "xy"
#    writeghostdata("../" + datafile)
#    Dir_seed = read_data("ghost_grain.data")
#    center = Get_Center_of_Mass(Dir_seed)
#    D, center = Get_Mass_Radius(Dir_seed, center, dimension, plane)
#    radius = math.ceil(D) + 5
#    print("Grain Radius:", D)  #'Coordinates cylinder: ', coord,
#    GD = D

#    nc = calc_polymer(mass, nm, polymer_type, p_fra, custom=mw)

# Get new cell for polymod
# Run polymod
#    shift, cell_polymod, cpos, caxis, clength = parameters(
#        periodic, cell_sizes, coord, boundaries, radius
#    )
#    input_polymod(nc, nm, cell_polymod, cpos, caxis, radius, clength, mctemp)

#    p_flag = 1
#    while p_flag:
#        print("Running PolymerModeler, generate initial configuration of polymers")
#        run_polymod()

#        shutil.copy("atoms.dat", "./str2lammps/types")
#        shutil.copy("bonds.dat", "./str2lammps/types")
#        p = polymer_type.split("-")[0]
#        shutil.copy(
#            "./str2lammps/types/%s_atom_type.dat" % p,
#            "./str2lammps/types/atom_type.dat",
#        )
#        shutil.copy("./str2lammps/types/%s_atom_type.dat" % p, "%s_atom_type.dat" % p)
#        shutil.copy("bond_type.dat", "./str2lammps/types")

#        addatomtype("%s_atom_type.dat" % p, "atom_type.dat", elementlist)

# Str2Lammps
#        os.chdir("./str2lammps/")
# print os.getcwd()
#        new_typing, data4lammps = run_data4lammps("Gasteiger", "Dreiding", cell_polymod)
# subprocess.Popen(new_typing, shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# print(data4lammps)
#        print("Generating initial polymer datafile")
#        os.system(data4lammps)
# return_code = subprocess.Popen(data4lammps, shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#        os.chdir("../unwrap/")

#        print("Unwrap polymers")
#        os.system("../lmp_mpi < uw.in")
# return_code = subprocess.Popen('../lmp_mpi < uw.in', shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#        dump2data("expanded.data", "unwrap.dump", "step0.data")
#        os.chdir("..")
#        os.system("cp unwrap/step0.data .")
#        os.system("cp unwrap/unwrap.dump .")
#        Dir_dump = read_dump("unwrap.dump")

#        na = 0
#        for t in Dir_dump:
#            na = Dir_dump[t]["NUMBER OF ATOMS"]
#            break

#        write_minimize(
#            "minimize.in",
#            "str2lammps/X6paircoeffs.txt",
#            "str2lammps/types/%s_atom_type.dat" % p,
#        )

#        if rnp_flag:
#            print("Calculating recommended number of processors")
#            ppn = grabprocessors("minimize.in")
#            print(ppn)
#            nodes, ppn, np = correctppn(ppn)
#            print("%d of processors will be in use" % np)

#        print("Running LAMMPS, minimize the intial configuration")
#        if np_flag:
#            os.system("mpiexec -np %d ./lmp_mpi < minimize.in" % np)
#        else:
#            os.system("./lmp_mpi < minimize.in")
# return_code = subprocess.Popen('./lmp_mpi < minimize.in', shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#        os.system("cp log.lammps log.minimize")

#        print("Runing 10 steps for coating")

#        center = [cpos[0], cpos[1], 0]

#        for i in range(10):
#            print("Step %d" % (i + 1))

#            input_coat4(
#                na, cpos, "step%d.data" % (i * 3 + 1), "step%d.data" % (i * 3 + 2)
#            )

#            print("Running LAMMPS, assign velocities")
#            os.system("./lmp_mpi < run_coat4.in")
# return_code = subprocess.Popen('./lmp_mpi < run_coat4.in', shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#            os.system("cp log.lammps log.run_coat4_%d" % i)

#            Dir_seed = read_data("step%d.data" % (i * 3 + 2))
#            oradius, center = Get_Mass_Radius(Dir_seed, center, dimension, plane)
#            print("Current outter radius:", oradius)

#            Dir = read_data("step%d.data" % (i * 3 + 2))
#            if "Pair Coeffs" in Dir:
#                del Dir["Pair Coeffs"]
#            write_data(Dir, "step%d.data" % (i * 3 + 3))

#            input_coat5(
#                cpos,
#                radius - 5,
#                oradius,
#                "step%d.data" % (i * 3 + 3),
#                "step%d.data" % (i * 3 + 4),
#                i * 5000,
#                "str2lammps/X6paircoeffs.txt",
#                mctemp,
#            )

#            print("Running LAMMPS, coat polymers")
#            if np_flag:
#                os.system("mpiexec -np %d ./lmp_mpi < run_coat5.in" % np)
#            else:
#                os.system("./lmp_mpi < run_coat5.in")
# return_code = subprocess.Popen('./lmp_mpi < run_coat5.in', shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#            os.system("cp log.lammps log.run_coat5_%d" % i)

#            if not os.path.exists("step%d.data" % (i * 3 + 4)):
#                break

#        if os.path.exists("step31.data"):
#            os.system("mv step31.data ../polymer.data")
#            p_flag = 0
#        else:
#            print("Failed. Generating another configuration")

#    print("Polymer datafile polymer.data is ready, now combining with grain datafile")

#    Dir_polymer = read_data("../polymer.data")
#    vector = [-cpos[0], -cpos[1], 0]
#    Dir_polymer = data_Translation(Dir_polymer, vector, box=1)
#    vector = [-coord["x"], -coord["y"], 0]
#    Dir_RDX = data_Translation(Dir_RDX, vector)
#    Dir_data = add_data(Dir_polymer, Dir_RDX)
#    write_data(Dir_data, "../initial_polymer_grain.data", v=0)
#    print("Initial combined datafile initial_polymer_grain.data is ready to use")
#    print("Center Coordinates: x, ", 0, "y, ", 0)
#    print("Calculating initial polymer thickness")
#    Dir_seed = read_data("../initial_polymer_grain.data")
#    center = [0, 0, 0]
#    D, center = Get_Mass_Radius(Dir_seed, center, dimension, plane)
#    PT = D - GD
#    print("Initial polymer thickness:", PT)
# des=open('p_thickness.dat','w')
# des.write(str(PT))
# des.close()
#    shutil.copy("atom_type.dat", "./str2lammps/types")
#    os.chdir("./str2lammps/")
#    new_typing, data4lammps = run_data4lammps("Gasteiger", "Dreiding", cell_polymod)
#    print("Generating new potentialfile")
#    os.system(data4lammps)
#    os.chdir("..")
#    shutil.copy("./str2lammps/X6paircoeffs.txt", "X6paircoeffs.txt")
#    shiftpotential("X6paircoeffs.txt", "../potential.mod", "potential.mod", Dir_polymer)
#    print("Runing 10 steps for coating with grain")

#    center = [0, 0, 0]

#    input_coat6(na, "../initial_polymer_grain.data", "step32.data")

#    if rnp_flag:
#        print("Calculating recommended number of processors")
#        ppn = grabprocessors("run_coat6.in")
#        print(ppn)
#        nodes, ppn, np = correctppn(ppn)
#        print("%d of processors will be in use" % np)

#    for i in range(10):
#        print("Step %d" % (i + 1))

#        if i:
#            input_coat6(na, "step%d.data" % (i * 3 + 31), "step%d.data" % (i * 3 + 32))

#        print("Running LAMMPS, assign velocities")
#        os.system("./lmp_mpi < run_coat6.in")
# return_code = subprocess.Popen('./lmp_mpi < run_coat4.in', shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#        os.system("cp log.lammps log.run_coat6_%d" % i)

#        Dir_seed = read_data("step%d.data" % (i * 3 + 32))
#        oradius, center = Get_Mass_Radius(Dir_seed, center, dimension, plane)
#        print("Current outter radius:", oradius)

#        input_coat7(
#            oradius,
#            "step%d.data" % (i * 3 + 32),
#            "step%d.data" % (i * 3 + 33),
#            i * 5000,
#            Dir_polymer["atom types"],
#            len(elementlist),
#            mctemp,
#        )

#        print("Running LAMMPS, coat polymers")
#        if np_flag:
#            os.system("mpiexec -np %d ./lmp_mpi < run_coat7.in" % np)
#        else:
#            os.system("./lmp_mpi < run_coat7.in")
# return_code = subprocess.Popen('./lmp_mpi < run_coat5.in', shell=True,
# stdin=subprocess.PIPE, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# stdout,stderr = return_code.communicate()

#        os.system("cp log.lammps log.run_coat7_%d" % i)

#        Dir = read_data("step%d.data" % (i * 3 + 33))
#        if "Pair Coeffs" in Dir:
#            del Dir["Pair Coeffs"]
#        if "Bond Coeffs" in Dir:
#            del Dir["Bond Coeffs"]
#        if "Angle Coeffs" in Dir:
#            del Dir["Angle Coeffs"]
#        if "Dihedral Coeffs" in Dir:
#            del Dir["Dihedral Coeffs"]
#        if "Improper Coeffs" in Dir:
#            del Dir["Improper Coeffs"]
#        write_data(Dir, "step%d.data" % (i * 3 + 34))

#    os.system("mv step61.data ../%s" % output)
#    print("Datafile %s is ready to use" % output)
#    print("Center Coordinates: x, ", 0, "y, ", 0)
#    print("Calculating initial polymer thickness")
#    Dir_seed = read_data("../%s" % output)
#    center = [0, 0, 0]
#    D, center = Get_Mass_Radius(Dir_seed, center, dimension, plane)
#    PT = D - GD
#    print("Polymer thickness:", PT)
#    des = open("p_thickness.dat", "w")
#    des.write(str(PT))
#    des.close()
#    return True


def zen(with_attribution=True):
    quote = """Beautiful is better than ugly.
    Explicit is better than implicit.
    Simple is better than complex.
    Complex is better than complicated.
    Flat is better than nested.
    Sparse is better than dense.
    Readability counts.
    Special cases aren't special enough to break the rules.
    Although practicality beats purity.
    Errors should never pass silently.
    Unless explicitly silenced.
    In the face of ambiguity, refuse the temptation to guess.
    There should be one-- and preferably only one --obvious way to do it.
    Although that way may not be obvious at first unless you're Dutch.
    Now is better than never.
    Although never is often better than *right* now.
    If the implementation is hard to explain, it's a bad idea.
    If the implementation is easy to explain, it may be a good idea.
    Namespaces are one honking great idea -- let's do more of those!"""

    if with_attribution:
        quote += "\n\tTim Peters"

    return quote


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
    infile = sys.argv[1]
    # args = polymerxtal.read_input(infile)
    # main(arg)
