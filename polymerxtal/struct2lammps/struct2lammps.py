"""
This module generates all necessary input files for LAMMPS simulations of molecular systems starting with an atomistic structure.
"""

import numpy as np
import os

from .createBonds import simulation_box, create_bonds
from .createDreidingTypes import create_dreiding_types
from .readFiles import convert_structure
from .runData4Lammps import run_data4lammps


def save_structure(file, outfilePDB=""):

    if not outfilePDB:
        if not os.path.exists("bonds"):
            os.mkdir("bonds")
        outfilePDB = "./bonds/pdbfile.pdb"
    structure_name = file
    os.rename(structure_name, structure_name.replace(" ", ""))
    structure_path = structure_name.replace(" ", "")
    # print('File uploaded succesfully')

    extension = os.path.splitext(structure_path)[1]
    if extension == ".cif":
        os.system(
            "obabel {0} -O {1} --fillUC strict --title convert ---errorlevel 1".format(
                structure_path, outfilePDB
            )
        )
    else:
        os.system(
            "obabel {0} -O {1} --title convert ---errorlevel 1".format(
                structure_path, outfilePDB
            )
        )

    # imolecule.draw("./bonds/pdbfile.pdb",size=(1000, 300))#, camera_type="orthographic")
    return outfilePDB


def conect_in_pdb(infile):

    with open(infile) as pdb_file:
        for line in pdb_file:
            if line.startswith("CONECT"):
                return True

        return False


def Create_Data_File(
    file,
    bondscale=1.1,
    lattice_in_file=False,
    ffield="Dreiding",
    charge="Gasteiger",
    xlo=0,
    xhi=50,
    ylo=0,
    yhi=50,
    zlo=0,
    zhi=50,
    alpha=90,
    beta=90,
    gamma=90,
):
    infile = save_structure(file)

    boundaries = np.zeros(9)
    boundaries[0] = xlo
    boundaries[1] = xhi
    boundaries[2] = ylo
    boundaries[3] = yhi
    boundaries[4] = zlo
    boundaries[5] = zhi
    boundaries[6] = alpha
    boundaries[7] = beta
    boundaries[8] = gamma

    extension = os.path.splitext(infile)[1]
    infile_is_pdb = False

    if extension == ".pdb" and conect_in_pdb(infile):
        print("Connectivity already in {}. Copying information ...".format(infile))
        os.system(f"cp {infile} ./bonds/connected_pdb.pdb")
        outfilelmpdat, outfilepdb = convert_structure(infile)
        boundaries = simulation_box(outfilelmpdat, lattice_in_file, boundaries)
        infile_is_pdb = True

    else:
        print("Creating bonds ...")
        outfilelmpdat, outfilepdb = convert_structure(infile)
        create_bonds(outfilelmpdat, outfilepdb, bondscale, lattice_in_file, boundaries)

    outputName = infile.split("/")[-1].split(".")[0]
    create_dreiding_types(
        "./bonds/connected_pdb.pdb", outputName, infile_is_pdb, ffield
    )
    print("Getting information...")
    run_data4lammps(charge, infile_is_pdb, boundaries, ffield)
    os.system("rm *index *name")
