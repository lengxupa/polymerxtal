#!/usr/bin/python
import os
import sys
from subprocess import call, Popen, PIPE
import glob
import string
import time
import re
import numpy as np
import datetime
import math


def copy_file(ifile, ofile, sentence, i):
    src = file(ifile, "r")
    des = file(ofile, "w")
    for line in src.readlines():
        if not (line.find(sentence)):
            line = sentence + i + "\n"
        des.write(line)
    src.close()
    des.close()


def change_box(ifile, ofile):
    find_box("new0.data")
    src = file(ifile, "r")
    des = file(ofile, "w")
    box = file("Box.txt", "r")
    for line in src.readlines():
        if line.find("xlo xhi") != -1:
            ln = box.readline()
            des.write("1 extra bond per atom\n")
            des.write("50 extra angle per atom\n")
            des.write("50 extra dihedral per atom\n")
            des.write("50 extra improper per atom\n")
            des.write("20 extra special per atom\n")
            des.write("\n")
            des.write(ln[0:-1] + " xlo xhi\n")
            continue
        elif line.find("ylo yhi") != -1:
            ln = box.readline()
            des.write(ln[0:-1] + " ylo yhi\n")
            continue
        elif line.find("zlo zhi") != -1:
            ln = box.readline()
            des.write(ln[0:-1] + " zlo zhi\n")
            continue
        des.write(line)
    box.close()
    src.close()
    des.close()


def change_kaka(ifile, ofile):
    find_kaka(ifile)
    src = file(ifile, "r")
    des = file(ofile, "w")
    box = file("Box.txt", "r")
    for line in src.readlines():
        if not (line.find("0.0 50 xlo xhi")):
            ln = box.readline()
            des.write(ln[0:-1] + " xlo xhi\n")
            continue
        elif not (line.find("0.0 50 ylo yhi")):
            ln = box.readline()
            des.write(ln[0:-1] + " ylo yhi\n")
            continue
        elif not (line.find("0.0 50 zlo zhi")):
            ln = box.readline()
            des.write(ln[0:-1] + " zlo zhi\n")
            continue
        des.write(line)
    box.close()
    src.close()
    des.close()


def find_box(ifile):
    src = file(ifile, "r")
    des = file("a.m", "w")
    a = 0
    des.write("A=[")
    for line in src.readlines():
        if a and line.find("Velocities"):
            des.write(line)
        if not (line.find("Atoms")):
            a = 1
        elif not (line.find("Velocities")):
            a = 0
    des.write("];")
    src.close()
    des.close()
    write_box()


def find_kaka(ifile):
    src = file(ifile, "r")
    des = file("a.m", "w")
    a = 0
    des.write("A=[")
    for line in src.readlines():
        if a and line.find("Bonds"):
            des.write(line)
        if not (line.find("Atoms")):
            a = 1
        elif not (line.find("Bonds")):
            a = 0
    des.write("];")
    src.close()
    des.close()
    os.system("matlab < kaka.m")


def write_box():
    src = file("box.m", "w")
    src.write("run a.m;\n")
    src.write("ne=size(A,1);\n")
    src.write("B=sortrows(A);\n")
    src.write("xmin=min(A(:,5))-25;\n")
    src.write("xmax=max(A(:,5))+25;\n")
    src.write("ymin=min(A(:,6))-25;\n")
    src.write("ymax=max(A(:,6))+25;\n")
    src.write("zmin=min(A(4,7),A(ne-5,7));\n")
    src.write("zmax=max(A(4,7),A(ne-5,7));\n")
    src.write("Box=[xmin xmax\n")
    src.write("ymin ymax\n")
    src.write("zmin zmax];\n")
    src.write("dlmwrite('Box.txt',Box,' ');\n")
    src.close()
    os.system("matlab < box.m")


def write_sublmp(ifile):
    src = file("sub.sh", "w")
    src.write("#!/bin/tcsh\n")
    src.write("#PBS -l walltime=0:10:00\n")
    src.write("#PBS -l nodes=1:ppn=16\n")
    src.write("#PBS -N before\n")
    a = findqueue(600, 16)
    while True:
        if (
            a == "debug"
            or a == "ncn"
            or a == "prism"
            or a == "standby"
            or a == "strachan"
        ):
            break
        a = findqueue(600, 16)
    src.write("#PBS -q %s\n" % a)
    src.write("\n")
    src.write("set echo\n")
    src.write("\n")
    src.write("cd $PBS_O_WORKDIR\n")
    src.write("module load intel impi/5.1.2.150 lammps/15Feb16\n")
    src.write("mpirun -np 4 lmp < %s \n" % ifile)
    src.write("\n")
    src.close()


def findqueue(i, n):
    if i <= 600:
        i = 600
    if i > 600 and i <= 1200:
        i = 1200
    if i > 1200 and i <= 1800:
        i = 1800
    if i > 1800 and i <= 3600:
        i = 3600
    if i / 3600.0 > 1 and i / 3600.0 <= 12:
        if int(i) % 3600 or int(i) != i:
            i = (int(i) / 3600 + 1) * 3600
        else:
            i = int(i)
    if i / 3600.0 > 12 and i / 3600.0 <= 24:
        i = 86400
    if i / 86400.0 > 1 and i / 86400.0 <= 15:
        if int(i) % 86400 or int(i) != i:
            i = (int(i) / 86400 + 1) * 86400
        else:
            i = int(i)
    if i / 86400.0 > 15:
        i = 1209600
    while True:
        cmd = "ssh conte /usr/pbs/bin/qlist"
        tst = Popen(cmd, shell=True, stdout=PIPE)
        outp, err = tst.communicate()
        print outp
        if not (err):
            break
    outplala = outp.split("\n")
    for line in outplala:
        ln = rdln(line)
        if ln:
            line = ln
        else:
            continue
        if ord(line[0]) >= 97 and ord(line[0]) <= 123:
            ln = line.replace(",", "")
            content = ln.split(" ", 1)
            raw = content[1].strip().lstrip().rstrip(",").split(" ", 1)
            raw2 = raw[1].strip().lstrip().rstrip(",").split(" ", 1)
            raw3 = raw2[1].strip().lstrip().rstrip(",").split(" ", 1)
            raw4 = raw3[1].strip().lstrip().rstrip(",").split(" ", 1)
            raw5 = raw4[1].strip().lstrip().rstrip(",").split(" ", 1)
            list1 = []
            if content[0] == "debug":
                debug = {
                    "Name": content[0],
                    "Total": raw[0],
                    "Queue": raw2[0],
                    "Run": raw3[0],
                    "Free": raw4[0],
                    "Max Walltime": walltimetoseconds(raw5[0]),
                }
            if content[0] == "ncn":
                ncn = {
                    "Name": content[0],
                    "Total": raw[0],
                    "Queue": raw2[0],
                    "Run": raw3[0],
                    "Free": raw4[0],
                    "Max Walltime": walltimetoseconds(raw5[0]),
                }
            if content[0] == "prism":
                prism = {
                    "Name": content[0],
                    "Total": raw[0],
                    "Queue": raw2[0],
                    "Run": raw3[0],
                    "Free": raw4[0],
                    "Max Walltime": walltimetoseconds(raw5[0]),
                }
            if content[0] == "standby":
                standby = {
                    "Name": content[0],
                    "Total": raw[0],
                    "Queue": raw2[0],
                    "Run": raw3[0],
                    "Free": raw4[0],
                    "Max Walltime": walltimetoseconds(raw5[0]),
                }
            if content[0] == "strachan":
                strachan = {
                    "Name": content[0],
                    "Total": raw[0],
                    "Queue": raw2[0],
                    "Run": raw3[0],
                    "Free": raw4[0],
                    "Max Walltime": walltimetoseconds(raw5[0]),
                }
    a = 1
    k = 0
    t = 0
    while a:
        k += 1
        for dirc in (debug, ncn, prism, standby, strachan):
            if i <= dirc["Max Walltime"]:
                if eval(dirc["Free"]) >= n and eval(dirc["Queue"]) == 0:
                    a = 0
                    if n != 16 and dirc == debug:
                        continue
                    return dirc["Name"]
                    break
                if (
                    k == 2
                    and eval(dirc["Free"]) >= n
                    and eval(dirc["Queue"]) < eval(dirc["Free"])
                ):
                    if dirc == debug and eval(dirc["Queue"]) + eval(dirc["Run"]) >= 64:
                        continue
                    a = 0
                    return dirc["Name"]
                    break
                if k == 3 and eval(dirc["Free"]) >= n:
                    if dirc == debug and eval(dirc["Queue"]) + eval(dirc["Run"]) >= 64:
                        continue
                    a = 0
                    return dirc["Name"]
                    break
                if k == 4 and eval(dirc["Queue"]) == 0:
                    a = 0
                    return dirc["Name"]
                    break
                if k == 5 and eval(dirc["Run"]) > t:
                    t = eval(dirc["Run"])
                if k == 6 and eval(dirc["Run"]) == t:
                    if dirc == debug and eval(dirc["Queue"]) == 64:
                        continue
                    a = 0
                    return dirc["Name"]
                    break


def rdln(i):
    content = i.strip().lstrip().rstrip(",").split("#", 1)
    return content[0]


def walltimetoseconds(i):
    content = i.split(":")
    if eval(content[0]) < 24:
        pt = datetime.datetime.strptime(i, "%H:%M:%S")
        total_seconds = pt.second + pt.minute * 60 + pt.hour * 3600
    else:
        total_seconds = (
            eval(content[2]) + eval(content[1]) * 60 + eval(content[0]) * 3600
        )
    return total_seconds


def distance(a, b):
    t = 0
    for i in range(3):
        t += (b[i] - a[i]) * (b[i] - a[i])
    return np.sqrt(t)


def writeTk(Atoms):
    des = open("T.txt", "w")
    a = vec(Atoms[4], Atoms[10])
    for i in range(198):
        b = vec(Atoms[i * 7 + 4], Atoms[i * 7 + 10])
        des.write("%d %f %f\n" % (i, k(a, b), distance(Atoms[4], Atoms[i * 7 + 4])))
    des.close()


def readdata(a):
    src = open(a, "r")
    A = 0
    Atoms = {}
    for line in src.readlines():
        ln = line.split()
        if ln:
            if ln[0] == "Atoms":
                A = 1
                continue
            if ln[0] == "Bonds":
                A = 0
                continue
            if A:
                Atoms[eval(ln[0])] = [eval(i) for i in ln[4:]]
    src.close()
    return Atoms


def valibond(a):
    Atoms = readdata(a)
    Bonds = readbond(a)
    for key in Bonds:
        if distance(Atoms[key[0]], Atoms[key[1]]) > 1.6:
            print "bond no"
            return False
            break
    print "bond yes"
    return True


def validate(a):
    Atoms = readdata(a)
    for i in range(1, len(Atoms)):
        for j in range(i + 1, len(Atoms) + 1):
            if distance(Atoms[i], Atoms[j]) < 1:
                print "no"
                return False
                break
    print "yes"
    return True


def k(a, b):
    return np.dot(a, b) / leng(a) / leng(b)


def vec(a, b):
    return [b[i] - a[i] for i in range(3)]


def leng(a):
    t = 0
    for i in range(3):
        t += a[i] * a[i]
    return np.sqrt(t)


def readbond(a):
    src = open(a, "r")
    A = 0
    Bonds = []
    for line in src.readlines():
        ln = line.split()
        if ln:
            if ln[0] == "Bonds":
                A = 1
                continue
            if ln[0] == "Angles":
                A = 0
                continue
            if A and eval(ln[1]) == 2:
                Bonds.append([eval(ln[2]), eval(ln[3])])
    src.close()
    return Bonds


des = open("table.txt", "w")
for tac in ["iso", "ata", "syn"]:
    des.write(tac + "\n")
    for i0 in range(1, 19) if tac == "iso" else range(1, 10):
        for j0 in range(
            (i0 if tac == "iso" else 10), (10 if tac == "iso" and i0 < 10 else 19)
        ):
            if tac == "ata":
                j0 = i0 + 12 if i0 < 7 else i0 + 3
            i = "m%d" % i0
            j = "m%d" % j0
            if not (valibond("../%s/%s/graphene.data" % (tac, i + j))):
                print tac, i, j
            if tac == "ata":
                break
            src = open("T%s.txt" % ("_".join([tac, i, j])))
            Dir = {}
            l = 1
            for line in src.readlines():
                ln = line.split()
                if eval(ln[0]) > 0 and eval(ln[1]) > 0.99:
                    Dir[eval(ln[0])] = eval(ln[2])
                if eval(ln[0]) > 0 and eval(ln[2]) < 1:
                    l = 0
                    break
            src.close()
            if l:
                for k in Dir:
                    if (k * 2 in Dir) and (k * 3 in Dir) and (k * 5 in Dir):
                        if -0.01 < (Dir[k * 2] - Dir[k] * 2) < 0.01:
                            if -0.01 < (Dir[k * 3] - Dir[k] * 3) < 0.01:
                                if -0.01 < (Dir[k * 5] - Dir[k] * 5) < 0.01:
                                    if validate(
                                        "../%s/%s/graphene.data" % (tac, i + j)
                                    ):
                                        print "here"
                                        des.write(i + " " + j + " " + ("%d\n" % k))
                                        break

des.close()
