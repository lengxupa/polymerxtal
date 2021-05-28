import os
import string

from .Gasteiger import getGasteiger_parameters, getGasteigerCharge
from .getForcefield import *  # noqa: F403
from .handleAtoms import Atomtypes, AtomsInfo, AtomLink
from .handleBonds import *  # noqa: F403
from .PCFF import (
    PCFF_getAngletypes,
    PCFF_getDihstypes,
    PCFF_getImpstypes,
    PCFF_getAtommass,
    PCFF_getPairCoeffs,
    PCFF_readPairCoeffs,
    PCFF_getBondCoeffs,
    PCFF_getAngleCoeffs,
    getBBCoeffs,
    getBACoeffs,
    PCFF_getDihsCoeffs,
    getMBTCoeffs,
    getEBTCoeffs,
    getATCoeffs,
    getAATCoeffs,
    getBB13Coeffs,
    PCFF_getImpsCoeffs,
    getAACoeffs,
)
from .qeq import Qeq_charge_equilibration


##############################################################################################
def checkAtomtype(inpfile):
    atomtypes = Atomtypes(inpfile)
    Forcefieldfile = getDreidingParamFile()  # noqa: F405
    typedefault = ("H_", "C_3", "N_3", "O_3", "F_", "S_3", "Cl", "I_", "Br_")

    atypes = []
    flag = 0

    fin = open(Forcefieldfile, "r")
    dataline = fin.readline()
    while dataline != "" and dataline != "\n" and flag == 0:
        words = dataline[0 : len(dataline) - 1]  # noqa: E203
        if str(words).upper() == "ATOMTYPES":
            flag = 1
            dataline = fin.readline()
            words = dataline[0 : len(dataline) - 1].split()  # noqa: E203
            while str(words[0]).upper() != "END":
                atype = str(words[0])
                atypes.append(atype)
                dataline = fin.readline()
                words = dataline[0 : len(dataline) - 1].split()  # noqa: E203
        dataline = fin.readline()
    fin.close()

    # print(atypes)
    anychange = "NO"
    changed = []
    for i in range(len(atomtypes)):
        atomtypeID = atomtypes[i][0]
        atomtype = atomtypes[i][1]
        if atomtype not in atypes:
            anychange = "YES"
            for j in range(len(typedefault)):
                deftype = typedefault[j]
                if atomtype[0:2] == deftype[0:2]:
                    atomtypes[i][1] = deftype
                    changed.append([atomtype, deftype])
    for i in range(len(atomtypes)):
        atomtypeID = atomtypes[i][0]
        atomtype = atomtypes[i][1]
        if atomtype not in atypes:
            anychange = "YES"
            for j in range(len(typedefault)):
                deftype = typedefault[j]
                if atomtype[0] == deftype[0]:
                    atomtypes[i][1] = deftype
                    changed.append([atomtype, deftype])
    for i in range(len(atomtypes)):
        atomtypeID = atomtypes[i][0]
        atomtype = atomtypes[i][1]
        if atomtype not in atypes:
            anychange = "YES"
            deftype = "C_3"
            atomtypes[i][1] = deftype
            changed.append([atomtype, deftype])
    if anychange == "YES":
        fout = open("atom_type_reassigned.dat", "w")
        for i in range(len(atomtypes)):
            atomtypeID = atomtypes[i][0]
            atomtype = atomtypes[i][1]
            fout.write(str(atomtypeID) + " " + atomtype + "\n")
            # print >> fout, atomtypeID, atomtype
        fout.close()
    if anychange == "YES":
        wout = open("Datafile_warnings1.txt", "w")
        wout.write(
            "##==============Warning: Force field parameters============================"
            + "\n"
        )
        # print >> wout, "##==============Warning: Force field parameters============================"
        wout.write("##  Atom type is re-assigned as following:" + "\n")
        # print >> wout, "##  Atom type is re-assigned as following:"
        wout.write("##" + str(changed) + "\n")
        # print >> wout, "##", changed
        wout.write(
            "##==============Warning: Force field parameters============================"
            + "\n"
        )
        # print >> wout, "##==============Warning: Force field parameters============================"
        wout.close()
    return atomtypes, changed


##############################################################################################
def printCoeffs(fout, ptitle, ptypes, ptypecoeffs):
    fout.write("\n")
    # print >>fout
    fout.write(ptitle + "\n")
    # print >>fout,ptitle
    fout.write("\n")
    # print >>fout
    for i in range(len(ptypes)):
        outline = ""
        for j in range(len(ptypecoeffs[i + 1])):
            outline = outline + str(ptypecoeffs[i + 1][j]) + " \t"
        fout.write(outline + "\n")
        # print >>fout,outline


##############################################################################################
# Write out force field parameters.
def outputDreidingCoeffs(fout, atomtypes, bondtypes, angletypes, dihstypes, impstypes):
    warning1 = getPairCoeffs(  # noqa: F405
        atomtypes
    )  # coeffs are in file "paircoeffs.txt"  # noqa: F405
    paircoeffs = readPairCoeffs()
    fout.write("\n")
    # print >>fout
    fout.write("Pair Coeffs" + "\n")
    # print >>fout,"Pair Coeffs"
    fout.write("\n")
    # print >>fout
    for i in range(len(paircoeffs)):
        fout.write(
            "%3i  %12.6f %12.6f   %s %s"
            % (
                paircoeffs[i][0],
                paircoeffs[i][1],
                paircoeffs[i][2],
                paircoeffs[i][3],
                paircoeffs[i][4],
            )
            + "\n"
        )
        # print >>fout,'%3i  %12.6f %12.6f   %s %s' % \
        #   (paircoeffs[i][0],paircoeffs[i][1],paircoeffs[i][2],paircoeffs[i][3],paircoeffs[i][4])
    fout.write("\n")
    # print >>fout
    fout.write("Bond Coeffs" + "\n")
    # print >>fout,"Bond Coeffs"
    fout.write("\n")
    # print >>fout
    bondcoeffs, warning2 = getBondCoeffs(bondtypes)  # noqa: F405
    for i in range(len(bondtypes)):
        fout.write(
            "%3i  %12.6f  %12.6f  %s%s%s%s"
            % (
                bondcoeffs[i][0],
                bondcoeffs[i][1],
                bondcoeffs[i][2],
                str("  # "),
                bondcoeffs[i][3],
                str(" "),
                bondcoeffs[i][4],
            )
            + "\n"
        )
        # print >>fout,'%3i  %12.6f  %12.6f  %s%s%s%s' % (bondcoeffs[i][0],bondcoeffs[i][1],bondcoeffs[i][2],str("  \
        #   # "),bondcoeffs[i][3],str(" "),bondcoeffs[i][4])
    fout.write("\n")
    # print >>fout
    fout.write("Angle Coeffs" + "\n")
    # print >>fout,"Angle Coeffs"
    fout.write("\n")
    # print >>fout
    anglecoeffs, warning3 = getAngleCoeffs(angletypes)  # noqa: F405
    for i in range(len(angletypes)):
        fout.write(
            "%3i  %12.6f  %12.6f  %s%s%s"
            % (
                anglecoeffs[i][0],
                anglecoeffs[i][1],
                anglecoeffs[i][2],
                str("  # X "),
                anglecoeffs[i][3],
                str(" X "),
            )
            + "\n"
        )
        # print >>fout,'%3i  %12.6f  %12.6f  %s%s%s' % (anglecoeffs[i][0],anglecoeffs[i][1],anglecoeffs[i][2],str("  #\
        #    X "),anglecoeffs[i][3],str(" X "))
    fout.write("\n")
    # print >>fout
    fout.write("Dihedral Coeffs" + "\n")
    # print >>fout,"Dihedral Coeffs"
    fout.write("\n")
    # print >>fout
    dihscoeffs, warning4 = getDihsCoeffs(dihstypes)  # noqa: F405
    for i in range(len(dihstypes)):
        fout.write(
            "%3i  %12.6f  %3i %3i %s%s%s%s%s%s%s%s"
            % (
                dihscoeffs[i][0],
                dihscoeffs[i][1],
                dihscoeffs[i][3],
                dihscoeffs[i][2],
                str("  # "),
                dihscoeffs[i][4],
                str(" "),
                dihscoeffs[i][5],
                str(" "),
                dihscoeffs[i][6],
                str(" "),
                dihscoeffs[i][7],
            )
            + "\n"
        )
        # print >>fout,'%3i  %12.6f  %3i %3i %s%s%s%s%s%s%s%s' % \
        #   (dihscoeffs[i][0],dihscoeffs[i][1],dihscoeffs[i][3],dihscoeffs[i][2],str("  # "),dihscoeffs[i][4],str(" \
        #   "),dihscoeffs[i][5],str(" "),dihscoeffs[i][6],str(" "),dihscoeffs[i][7])
    fout.write("\n")
    # print >>fout
    fout.write("Improper Coeffs" + "\n")
    # print >>fout,"Improper Coeffs"
    fout.write("\n")
    # print >>fout
    impscoeffs, warning5 = getImpsCoeffs(impstypes)  # noqa: F405
    for i in range(len(impscoeffs)):
        fout.write(
            "%3i  %12.6f  %12.6f  %s%s%s"
            % (
                impscoeffs[i][0],
                impscoeffs[i][1],
                impscoeffs[i][2],
                str("  #  "),
                impscoeffs[i][3],
                str(" X X X "),
            )
            + "\n"
        )
        # print >>fout,'%3i  %12.6f  %12.6f  %s%s%s' % (impscoeffs[i][0],impscoeffs[i][1],impscoeffs[i][2],str("  #  \
        #   "),impscoeffs[i][3],str(" X X X "))

    if (
        warning1 != ""
        or warning2 != ""
        or warning3 != ""
        or warning4 != ""
        or warning5 != ""
    ):
        wout = open("Datafile_warnings2.txt", "w")
        wout.write("##" + warning1 + "\n")
        # print >>wout,"##",warning1
        wout.write("##" + warning2 + "\n")
        # print >>wout,"##",warning2
        wout.write("##" + warning3 + "\n")
        # print >>wout,"##",warning3
        wout.write("##" + warning4 + "\n")
        # print >>wout,"##",warning4
        wout.write("##" + warning5 + "\n")
        # print >>wout,"##",warning5
        wout.write(
            "##==============Warning: Force field parameters============================"
            + "\n"
        )
        # print >>wout,"##==============Warning: Force field parameters============================"
        wout.close()


##############################################################################################
# Write out force field parameters.
def outputPCFFCoeffs(fout, atomtypes, bondtypes, angletypes, dihstypes, impstypes):
    PCFF_getPairCoeffs(atomtypes)  # coeffs are in file "paircoeffs.txt"
    paircoeffs = PCFF_readPairCoeffs()
    printCoeffs(fout, "Pair Coeffs", atomtypes, paircoeffs)

    bondcoeffs = PCFF_getBondCoeffs(bondtypes)
    printCoeffs(fout, "Bond Coeffs", bondtypes, bondcoeffs)

    anglecoeffs = PCFF_getAngleCoeffs(angletypes)
    printCoeffs(fout, "Angle Coeffs", angletypes, anglecoeffs)

    BBcoeffs = getBBCoeffs(angletypes, bondtypes, bondcoeffs)
    printCoeffs(fout, "BondBond Coeffs", angletypes, BBcoeffs)

    BAcoeffs = getBACoeffs(angletypes, bondtypes, bondcoeffs)
    printCoeffs(fout, "BondAngle Coeffs", angletypes, BAcoeffs)

    dihscoeffs = PCFF_getDihsCoeffs(dihstypes)
    printCoeffs(fout, "Dihedral Coeffs", dihstypes, dihscoeffs)

    MBTcoeffs = getMBTCoeffs(dihstypes, bondtypes, bondcoeffs)
    printCoeffs(fout, "MiddleBondTorsion Coeffs", dihstypes, MBTcoeffs)

    EBTcoeffs = getEBTCoeffs(dihstypes, bondtypes, bondcoeffs)
    printCoeffs(fout, "EndBondTorsion Coeffs", dihstypes, EBTcoeffs)

    ATcoeffs = getATCoeffs(dihstypes, angletypes, anglecoeffs)
    printCoeffs(fout, "AngleTorsion Coeffs", dihstypes, ATcoeffs)

    AATcoeffs = getAATCoeffs(dihstypes, angletypes, anglecoeffs)
    printCoeffs(fout, "AngleAngleTorsion Coeffs", dihstypes, AATcoeffs)

    BB13coeffs = getBB13Coeffs(dihstypes, bondtypes, bondcoeffs)
    printCoeffs(fout, "BondBond13 Coeffs", dihstypes, BB13coeffs)

    impscoeffs = PCFF_getImpsCoeffs(impstypes)
    printCoeffs(fout, "Improper Coeffs", impstypes, impscoeffs)

    AAcoeffs = getAACoeffs(impstypes, angletypes, anglecoeffs)
    printCoeffs(fout, "AngleAngle Coeffs", impstypes, AAcoeffs)


########################################################################
def getFileName():
    fin = open("structure.name", "r")
    dataline = fin.readline()
    words = dataline[0 : len(dataline) - 1].split()  # noqa: E203
    structureName = words[0]
    return structureName


########################################################################
def getForcefield():
    fin = open("forcefield.name", "r")
    dataline = fin.readline()
    words = dataline[0 : len(dataline) - 1].split()  # noqa: E203
    forcefieldName = words[0]
    return forcefieldName


########################################################################
def readPairCoeffs():
    paircoeffs = []

    fin = open("LJpaircoeffs.txt", "r")
    dataline = fin.readline()
    while dataline != "":
        words = dataline[0 : len(dataline) - 1].split()  # noqa: E203
        atomtype = eval(words[1])
        D0 = eval(words[2])
        R0 = eval(words[3])
        C1 = words[4]
        C2 = words[5]
        paircoeffs.append([atomtype, D0, R0, C1, C2])
        dataline = fin.readline()
    return paircoeffs


########################################################################
# Assumes cubic cell
def createReaxDatafile(
    forcefield, structureName, xlo, xhi, ylo, yhi, zlo, zhi, chargeMethod
):
    inpfile = "atom_type.dat"
    if str(forcefield).upper() == "DREIDING":
        atomtypes, anychange = checkAtomtype(inpfile)
    else:
        atomtypes = Atomtypes(inpfile)
        # anychange = []

    print("Atomtypes total=", len(atomtypes))
    inpfile = "atoms.dat"
    baseatoms = AtomsInfo(inpfile)
    print("Atoms total=", len(baseatoms))
    natomtype = len(atomtypes)
    totalatoms = len(baseatoms)

    atommass = getAtommass(atomtypes)  # noqa: F405
    ####################################################################
    # Output reaxFF data file for lammps
    datafile = structureName + "_reaxFF.data"
    fout = open(datafile, "w")
    fout.write("LAMMPS data file using " + forcefield + " for " + structureName + "\n")
    # print >>fout,"LAMMPS data file using "+forcefield+" for "+structureName
    fout.write("\n")
    # print >>fout
    fout.write(str(totalatoms) + "  atoms" + "\n")
    # print >>fout,str(totalatoms)+"  atoms"
    fout.write(str(natomtype) + "  atom types" + "\n")
    # print >>fout,str(natomtype)+"  atom types"
    fout.write("\n")
    # print >>fout
    fout.write(xlo + " " + xhi + " xlo xhi" + "\n")
    # print >>fout,xlo+" "+xhi+" xlo xhi"
    fout.write(ylo + " " + yhi + " ylo yhi" + "\n")
    # print >>fout,ylo+" "+yhi+" ylo yhi"
    fout.write(zlo + " " + zhi + " zlo zhi" + "\n")
    # print >>fout,zlo+" "+zhi+" zlo zhi"
    fout.write("\n")
    # print >>fout
    fout.write("Masses" + "\n")
    # print >>fout,"Masses"
    fout.write("\n")
    # print >>fout
    for i in range(len(atommass)):
        atomtype = atommass[i][2]
        # atomtype=atomtype[0]
        fout.write(
            "%3i  %12.6f   %s%s"
            % (atommass[i][0], atommass[i][1], str("  # "), atomtype)
            + "\n"
        )
        # print >>fout,'%3i  %12.6f   %s%s' % (atommass[i][0],atommass[i][1],str("  # "), atomtype)

    ####################################################################
    # Output atom data
    fout.write("\n")
    # print >>fout
    fout.write("Atoms # full" + "\n")
    # print >>fout, "Atoms"
    fout.write("\n")
    # print >>fout
    for i in range(len(baseatoms)):
        dataline = str(baseatoms[i])
        w = string.split(dataline[1 : len(dataline) - 1], ",")  # noqa: E203
        fout.write(
            (
                "%6d %3d %3d %10.5f %15.8f %15.8f %15.8f"
                % (
                    eval(w[0]),
                    eval(w[1]),
                    eval(w[2]),
                    eval(w[3]),
                    eval(w[4]),
                    eval(w[5]),
                    eval(w[6]),
                )
            )
            + "\n"
        )
        # print >>fout, ('%6d %3d %3d %10.5f %15.8f %15.8f %15.8f' %
        #               (eval(w[0]),eval(w[1]),eval(w[2]),eval(w[3]),eval(w[4]),eval(w[5]),eval(w[6])))
    fout.write("\n")
    # print >>fout
    fout.close()
    print(datafile + " created!")
    return datafile
    ####################################################################


def createDatafile(
    forcefield, structureName, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, chargeMethod
):
    inpfile = ".tmp/types/newatom_type.dat"
    if str(forcefield).upper() == "DREIDING":
        atomtypes, anychange = checkAtomtype(inpfile)
    else:
        atomtypes = Atomtypes(inpfile)
        anychange = []
    inpfile = ".tmp/types/newatoms.dat"
    baseatoms = AtomsInfo(inpfile)
    # Update bondtype if default types used
    inpfile = ".tmp/types/newbond_type.dat"
    bondtypes = getBondtypes(inpfile)  # noqa: F405
    for i in range(len(bondtypes)):
        atom1type = bondtypes[i][1]
        atom2type = bondtypes[i][2]
        for j in range(len(anychange)):
            replaced = anychange[j][0]
            defatype = anychange[j][1]
            if atom1type.upper() == replaced.upper():
                bondtypes[i][1] = defatype
            if atom2type.upper() == replaced.upper():
                bondtypes[i][2] = defatype
    inpfile = ".tmp/types/newbonds.dat"
    basebonds = getBonds(inpfile, 0, 1)  # noqa: F405

    print("Equilibrating charge... \n")

    if chargeMethod == "Gasteiger":
        # Replace charge for Gasteiger charge
        forcefield = getForcefield()
        Gparas = getGasteiger_parameters(forcefield)
        Q = getGasteigerCharge(Gparas, atomtypes, baseatoms, basebonds)

    elif chargeMethod == "QEq":
        Q = Qeq_charge_equilibration(baseatoms)

    atomlinks = AtomLink(baseatoms, basebonds)
    # print("Links generated")

    baseangles = createAngles(atomlinks)  # noqa: F405
    # print("Angles generated")
    if str(forcefield).upper() == "PCFF":
        angletypes, baseangles = PCFF_getAngletypes(baseangles, baseatoms, atomtypes)
    if str(forcefield).upper() == "DREIDING":
        angletypes, baseangles = getAngletypes(  # noqa: F405
            baseangles, baseatoms, atomtypes
        )  # noqa: F405
    # print("Angles updated")

    basedihs = createDihedrals(atomlinks, basebonds)  # noqa: F405
    # print("Dihs generated")
    if str(forcefield).upper() == "PCFF":
        dihstypes, basedihs = PCFF_getDihstypes(basedihs, baseatoms, atomtypes)
    if str(forcefield).upper() == "DREIDING":
        dihstypes, basedihs = getDihstypes(basedihs, baseatoms, atomtypes)  # noqa: F405
    # print("Dihs updated")

    baseimps = createImpropers(atomlinks)  # noqa: F405
    # print("Imps generated")
    if str(forcefield).upper() == "PCFF":
        impstypes, baseimps = PCFF_getImpstypes(baseimps, baseatoms, atomtypes)
    if str(forcefield).upper() == "DREIDING":
        impstypes, baseimps = getImpstypes(baseimps, baseatoms, atomtypes)  # noqa: F405
    # print("Imps updated")

    ####################################################################
    # Total quantities
    natomtype = len(atomtypes)
    nbondtype = len(bondtypes)
    nangletype = len(angletypes)
    ndihstype = len(dihstypes)
    nimpstype = len(impstypes)

    totalatoms = len(baseatoms)
    totalbonds = len(basebonds)
    totalangles = len(baseangles)
    totaldihs = len(basedihs)
    totalimps = len(baseimps)
    ####################################################################
    atommass = getAtommass(atomtypes)  # noqa: F405
    if str(forcefield).upper() == "PCFF":
        atommass = PCFF_getAtommass(atomtypes)

    # Output Lammps data file
    ####################################################################
    # Output head of data file for lammps
    datafile = structureName + ".data"
    # datafile = "LAMMPSDataFile.data"
    fout = open(datafile, "w")
    fout.write("LAMMPS data file using " + forcefield + " for " + structureName + "\n")
    # print >>fout,"LAMMPS data file using "+forcefield+" for "+structureName
    fout.write("\n")
    # print >>fout
    fout.write(str(totalatoms) + "  atoms" + "\n")
    # print >>fout,str(totalatoms)+"  atoms"
    fout.write(str(totalbonds) + "  bonds" + "\n")
    # print >>fout,str(totalbonds)+"  bonds"
    fout.write(str(totalangles) + "  angles" + "\n")
    # print >>fout,str(totalangles)+"  angles"
    fout.write(str(totaldihs) + "  dihedrals" + "\n")
    # print >>fout,str(totaldihs)+"  dihedrals"
    fout.write(str(totalimps) + "  impropers" + "\n")
    # print >>fout,str(totalimps)+"  impropers"
    fout.write("\n")
    # print >>fout
    fout.write(str(natomtype) + "  atom types" + "\n")
    # print >>fout,str(natomtype)+"  atom types"
    fout.write(str(nbondtype) + "  bond types" + "\n")
    # print >>fout,str(nbondtype)+"  bond types"
    fout.write(str(nangletype) + "  angle types" + "\n")
    # print >>fout,str(nangletype)+"  angle types"
    fout.write(str(ndihstype) + "  dihedral types" + "\n")
    # print >>fout,str(ndihstype)+"  dihedral types"
    fout.write(str(nimpstype) + "  improper types" + "\n")
    # print >>fout,str(nimpstype)+"  improper types"
    fout.write("\n")
    # print >>fout
    fout.write(xlo + " " + xhi + " xlo xhi" + "\n")
    # print >>fout,xlo+" "+xhi+" xlo xhi"
    fout.write(ylo + " " + yhi + " ylo yhi" + "\n")
    # print >>fout,ylo+" "+yhi+" ylo yhi"
    fout.write(zlo + " " + zhi + " zlo zhi" + "\n")
    # print >>fout,zlo+" "+zhi+" zlo zhi"
    if xy == "0.0" and xz == "0.0" and yz == "0.0":
        fout.write("\n")
        # print >>fout
    else:
        fout.write(xy + " " + xz + " " + yz + " xy xz yz" + "\n")
        # print >>fout,xy+" "+xz+" "+yz+" xy xz yz"
    fout.write("\n")
    # print >>fout
    fout.write("Masses" + "\n")
    # print >>fout,"Masses"
    fout.write("\n")
    # print >>fout
    for i in range(len(atommass)):
        fout.write(
            "%3i  %12.6f   %s%s"
            % (atommass[i][0], atommass[i][1], str("  # "), atommass[i][2])
            + "\n"
        )
        # print >>fout,'%3i  %12.6f   %s%s' % (atommass[i][0],atommass[i][1],str("  # "), atommass[i][2])

    ####################################################################
    # Output data
    fout.write("\n")
    # print >>fout
    fout.write("Atoms # full" + "\n")
    # print >>fout, "Atoms"
    fout.write("\n")
    # print >>fout
    for i in range(len(baseatoms)):
        dataline = str(baseatoms[i])
        w = dataline[1 : len(dataline) - 1].split(",")  # noqa: E203
        fout.write(
            (
                "%6d %3d %3d %10.5f %15.8f %15.8f %15.8f %3d %3d %3d"
                % (
                    eval(w[0]),
                    eval(w[1]),
                    eval(w[2]),
                    Q[i + 1],
                    eval(w[4]),
                    eval(w[5]),
                    eval(w[6]),
                    eval(w[7]),
                    eval(w[8]),
                    eval(w[9]),
                )
            )
            + "\n"
        )
        # print >>fout, ('%6d %3d %3d %10.5f %15.8f %15.8f %15.8f %3d %3d %3d' %
        #               (eval(w[0]),eval(w[1]),eval(w[2]),Q[i+1],eval(w[4]),eval(w[5]),eval(w[6]), eval(w[7]), \
        #                   eval(w[8]), eval(w[9])))
    fout.write("\n")
    # print >>fout
    fout.write("Bonds" + "\n")
    # print >>fout, "Bonds"
    fout.write("\n")
    # print >>fout
    for i in range(len(basebonds)):
        dataline = str(basebonds[i])
        words = dataline[1 : len(dataline) - 1].split(",")  # noqa: E203
        outline = ""
        for i in range(len(words)):
            outline = outline + str(words[i]) + " "
        fout.write(outline + "\n")
        # print >>fout, outline
    fout.write("\n")
    # print >>fout
    fout.write("Angles" + "\n")
    # print >>fout, "Angles"
    fout.write("\n")
    # print >>fout
    for i in range(len(baseangles)):
        dataline = str(baseangles[i])
        words = dataline[1 : len(dataline) - 1].split(",")  # noqa: E203
        outline = ""
        for i in range(len(words)):
            outline = outline + str(words[i]) + " "
        fout.write(outline + "\n")
        # print >>fout, outline
    fout.write("\n")
    # print >>fout
    fout.write("Dihedrals" + "\n")
    # print >>fout, "Dihedrals"
    fout.write("\n")
    # print >>fout
    for i in range(len(basedihs)):
        dataline = str(basedihs[i])
        words = dataline[1 : len(dataline) - 1].split(",")  # noqa: E203
        outline = ""
        for i in range(len(words)):
            outline = outline + str(words[i]) + " "
        fout.write(outline + "\n")
        # print >>fout, outline
    fout.write("\n")
    # print >>fout
    fout.write("Impropers" + "\n")
    # print >>fout, "Impropers"
    fout.write("\n")
    # print >>fout
    for i in range(len(baseimps)):
        dataline = str(baseimps[i])
        words = dataline[1 : len(dataline) - 1].split(",")  # noqa: E203
        outline = ""
        for i in range(len(words)):
            outline = outline + str(words[i]) + " "
        fout.write(outline + "\n")
        # print >>fout, outline
    # Coeffs
    if str(forcefield).upper() == "DREIDING":
        outputDreidingCoeffs(
            fout, atomtypes, bondtypes, angletypes, dihstypes, impstypes
        )
    if str(forcefield).upper() == "PCFF":
        outputPCFFCoeffs(fout, atomtypes, bondtypes, angletypes, dihstypes, impstypes)
    fout.close()

    print(datafile + " created!")
    if os.path.exists("Datafile_warnings.txt"):
        cmd2 = "rm Datafile_warnings.txt"
        os.system(cmd2)
    if os.path.exists("Datafile_warnings1.txt"):
        cmd1 = "cat Datafile_warnings1.txt >>Datafile_warnings.txt"
        cmd2 = "rm Datafile_warnings1.txt"
        os.system(cmd1)
        os.system(cmd2)
    if os.path.exists("Datafile_warnings2.txt"):
        cmd1 = "cat Datafile_warnings2.txt >>Datafile_warnings.txt"
        cmd2 = "rm Datafile_warnings2.txt"
        os.system(cmd1)
        os.system(cmd2)
    return datafile


########################################################################
