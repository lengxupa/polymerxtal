"""
buttondef.py
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer, tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline polymer simulations.

Handles functions to build molecular-level polymer crystal structures
"""

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
    src = file(ifile, 'r')
    des = file(ofile, 'w')
    for line in src.readlines():
        if (not (line.find(sentence))):
            line = sentence + i + '\n'
        des.write(line)
    src.close()
    des.close()


def change_box(ifile, ofile):
    find_box('new0.data')
    src = file(ifile, 'r')
    des = file(ofile, 'w')
    box = file('Box.txt', 'r')
    for line in src.readlines():
        if (line.find('xlo xhi') != -1):
            ln = box.readline()
            des.write('1 extra bond per atom\n')
            des.write('50 extra angle per atom\n')
            des.write('50 extra dihedral per atom\n')
            des.write('50 extra improper per atom\n')
            des.write('20 extra special per atom\n')
            des.write('\n')
            des.write(ln[0:-1] + ' xlo xhi\n')
            continue
        elif (line.find('ylo yhi') != -1):
            ln = box.readline()
            des.write(ln[0:-1] + ' ylo yhi\n')
            continue
        elif (line.find('zlo zhi') != -1):
            ln = box.readline()
            des.write(ln[0:-1] + ' zlo zhi\n')
            continue
        des.write(line)
    box.close()
    src.close()
    des.close()


def change_kaka(ifile, ofile):
    find_kaka(ifile)
    src = file(ifile, 'r')
    des = file(ofile, 'w')
    box = file('Box.txt', 'r')
    for line in src.readlines():
        if (not (line.find('0.0 50 xlo xhi'))):
            ln = box.readline()
            des.write(ln[0:-1] + ' xlo xhi\n')
            continue
        elif (not (line.find('0.0 50 ylo yhi'))):
            ln = box.readline()
            des.write(ln[0:-1] + ' ylo yhi\n')
            continue
        elif (not (line.find('0.0 50 zlo zhi'))):
            ln = box.readline()
            des.write(ln[0:-1] + ' zlo zhi\n')
            continue
        des.write(line)
    box.close()
    src.close()
    des.close()


def find_box(ifile):
    src = file(ifile, 'r')
    des = file('a.m', 'w')
    a = 0
    des.write('A=[')
    for line in src.readlines():
        if (a and line.find('Velocities')):
            des.write(line)
        if (not (line.find('Atoms'))):
            a = 1
        elif (not (line.find('Velocities'))):
            a = 0
    des.write('];')
    src.close()
    des.close()
    write_box()


def find_kaka(ifile):
    src = file(ifile, 'r')
    des = file('a.m', 'w')
    a = 0
    des.write('A=[')
    for line in src.readlines():
        if (a and line.find('Bonds')):
            des.write(line)
        if (not (line.find('Atoms'))):
            a = 1
        elif (not (line.find('Bonds'))):
            a = 0
    des.write('];')
    src.close()
    des.close()
    os.system('matlab < kaka.m')


def write_box():
    src = file('box.m', 'w')
    src.write('run a.m;\n')
    src.write('ne=size(A,1);\n')
    src.write('B=sortrows(A);\n')
    src.write('xmin=min(A(:,5))-25;\n')
    src.write('xmax=max(A(:,5))+25;\n')
    src.write('ymin=min(A(:,6))-25;\n')
    src.write('ymax=max(A(:,6))+25;\n')
    src.write('zmin=min(B(1,7),B(ne-5,7));\n')
    src.write('zmax=max(B(1,7),B(ne-5,7));\n')
    src.write('Box=[xmin xmax\n')
    src.write('ymin ymax\n')
    src.write('zmin zmax];\n')
    src.write("dlmwrite('Box.txt',Box,' ');\n")
    src.close()
    os.system('matlab < box.m')


def write_sublmp(seconds, N, ifile, tmp=0, te=''):
    src = file('sub.sh', "w")
    src.write("#!/bin/sh -l\n")
    src.write("# FILENAME:  sub.sh\n")
    walltime = secondstowalltime(seconds)
    ppn = grabprocessors(ifile)
    nodes, ppn, np = correctppn(ppn)
    if np > 100:
        a = 'standby'
        cluster = 'halstead'
        if seconds > 14400:
            walltime = '4:00:00'
    else:
        a, cluster = findqueue(seconds, np)
        while True:
            if (a == 'debug' or a == 'ncn' or a == 'prism' or a == 'standby' or a == 'strachan'):
                break
            a, cluster = findqueue(seconds, np)
        if a == 'standby' and seconds > 14400:
            walltime = '4:00:00'
        if tmp and seconds > 14400:
            start = time.time()
            stop = time.time()
            while (stop - start) < 75:
                if (a == 'ncn' or a == 'prism' or a == 'strachan'):
                    break
                a, cluster = findqueue(seconds, np)
                stop = time.time()
            walltime = secondstowalltime(seconds)
            if not (a == 'ncn' or a == 'prism' or a == 'strachan'):
                a = 'ncn'
    src.write("#SBATCH -t %s\n" % walltime)
    src.write("#SBATCH -N %d\n" % nodes)
    src.write("#SBATCH -n %d\n" % ppn)
    src.write("#SBATCH --job-name=%s\n" % N)
    src.write("#SBATCH -A %s\n" % a)
    src.write("\n")
    src.write("set echo\n")
    src.write("\n")
    src.write("cd $SLURM_SUBMIT_DIR\n")
    if cluster != 'conte':
        src.write("module load intel/17.0.1.132 impi/2017.1.132 lammps/7Aug19\n")
    else:
        src.write("module load intel impi/5.1.2.150 lammps/15Feb16\n")
    src.write("mpiexec -np %d lmp < %s\n" % (np, ifile))
    src.write("module load anaconda/5.3.1-py27\n")
    src.write("python nptsub.py\n")
    src.write("\n")
    src.close()
    return cluster


def correctppn(ppn):
    nodes = ppn / 20 + 1 if (ppn / 20 and ppn % 20) else ppn / 20 if ppn / 20 else 1
    if not ppn / 20:
        ppn = 20 if ppn > 16 else 16 if ppn > 10 else 10 if ppn > 8 else 8 if ppn > 4 else 4 if ppn > 2 else 2 if ppn > 1 else 1
    return nodes, 20 if ppn / 20 else ppn, nodes * 20 if ppn / 20 else ppn


def grabprocessors(ifile):
    src = open(ifile)
    for line in src.readlines():
        ln = line.split()
        if ln and ln[0] == 'read_data':
            datafile = ln[1]
            break
    src.close()
    Dir = read_data(datafile)
    atoms = Dir['atoms']
    if 'bonds' in Dir:
        return atoms / 2000 + 1
    else:
        return atoms / 1000 + 1


def read_data(ifile):
    Dir = {}
    Box = {}
    Masses = {}
    src = open(ifile, 'r')
    a = 1
    s = ''
    for line in src.readlines():
        if a:
            a = 0
            continue
        ln = line.split('#')[0].split()
        if ln:
            if len(ln) > 1 and ln[0].isdigit() and (not ln[1].isdigit()) and (not s):
                Dir[' '.join(ln[1:])] = eval(ln[0])
            if len(ln) == 4 and ln[2][1:] == 'lo' and ln[3][1:] == 'hi':
                Dir[ln[2]] = eval(ln[0])
                Dir[ln[3]] = eval(ln[1])
            if not (ln[0][0].isdigit() or ln[0][0] == '-'):
                s = ' '.join(ln)
                Dir[s] = {}
            if s and (ln[0][0].isdigit() or ln[0][0] == '-'):
                Dir[s][eval(ln[0])] = [eval(i) for i in ln[1:]]
    src.close()
    return Dir


def secondstowalltime(i):
    i = correction(i)
    return ("%d:%d0:00" % (i / 3600, i % 3600 / 600))


def correction(i):
    if (i <= 600):
        i = 600
    if (i > 600 and i <= 1200):
        i = 1200
    if (i > 1200 and i <= 1800):
        i = 1800
    if (i > 1800 and i <= 3600):
        i = 3600
    if (i / 3600.0 > 1 and i / 3600.0 <= 12):
        if (int(i) % 3600 or int(i) != i):
            i = (int(i) / 3600 + 1) * 3600
        else:
            i = int(i)
    if (i / 3600.0 > 12 and i / 3600.0 <= 24):
        i = 86400
    if (i / 86400.0 > 1 and i / 86400.0 <= 14):
        if (int(i) % 86400 or int(i) != i):
            i = (int(i) / 86400 + 1) * 86400
        else:
            i = int(i)
    if (i / 86400.0 > 14):
        i = 1209600
    return i


def findqueue(i, n):
    i = correction(i)
    Dir = {}
    for cluster in ['halstead', 'rice']:
        Dir[cluster] = {}
        while True:
            cmd = ("ssh %s /usr/bin/slist" % cluster)
            tst = Popen(cmd, shell=True, stdout=PIPE)
            outp, err = tst.communicate()
            print outp
            if (not (err)):
                break
        outplala = outp.split('\n')
        for line in outplala:
            ln = rdln(line)
            if (ln):
                line = ln
            else:
                continue
            if (ord(line[0]) >= 97 and ord(line[0]) <= 123):
                ln = line.replace(',', '')
                content = ln.split(' ', 1)
                raw = content[1].strip().lstrip().rstrip(',').split(' ', 1)
                raw2 = raw[1].strip().lstrip().rstrip(',').split(' ', 1)
                raw3 = raw2[1].strip().lstrip().rstrip(',').split(' ', 1)
                raw4 = raw3[1].strip().lstrip().rstrip(',').split(' ', 1)
                raw5 = raw4[1].strip().lstrip().rstrip(',').split(' ', 1)
                list1 = []
                if content[0] == 'debug' or content[0] == 'ncn' or content[0] == 'standby' or content[0] == 'strachan':
                    Dir[cluster][content[0]] = {
                        'Name': content[0],
                        'Total': raw[0],
                        'Queue': raw2[0],
                        'Run': raw3[0],
                        'Free': raw4[0],
                        'Max Walltime': walltimetoseconds(raw5[0])
                    }
                    if cluster == 'conte':
                        Dir[cluster][content[0]]['Max Walltime'] = 24 * 60 * 60
    a = 1
    k = 0
    t = 0
    while (a):
        k += 1
        for cluster in ['halstead', 'rice']:
            for dirc in Dir[cluster]:
                if (i <= Dir[cluster][dirc]['Max Walltime']):
                    if (eval(Dir[cluster][dirc]['Free']) >= n and eval(Dir[cluster][dirc]['Queue']) == 0):
                        if (n != 16 and dirc == 'debug'):
                            continue
                        a = 0
                        return Dir[cluster][dirc]['Name'], cluster
                        break
                    if (k == 2 and eval(Dir[cluster][dirc]['Free']) >= n
                            and eval(Dir[cluster][dirc]['Queue']) < eval(Dir[cluster][dirc]['Free'])):
                        if (n != 16 and dirc == 'debug'):
                            continue
                        if (dirc == 'debug'
                                and eval(Dir[cluster][dirc]['Queue']) + eval(Dir[cluster][dirc]['Run']) >= 64):
                            continue
                        a = 0
                        return Dir[cluster][dirc]['Name'], cluster
                        break
                    if (k == 3 and eval(Dir[cluster][dirc]['Free']) >= n):
                        if (n != 16 and dirc == 'debug'):
                            continue
                        if (dirc == 'debug'
                                and eval(Dir[cluster][dirc]['Queue']) + eval(Dir[cluster][dirc]['Run']) >= 64):
                            continue
                        a = 0
                        return Dir[cluster][dirc]['Name'], cluster
                        break
                    if (k == 4 and eval(Dir[cluster][dirc]['Queue']) == 0):
                        if (n != 16 and dirc == 'debug'):
                            continue
                        a = 0
                        return Dir[cluster][dirc]['Name'], cluster
                        break
                    if k == 5:
                        a = 0
                        break
    return 'standby', 'halstead'


def rdln(i):
    content = i.strip().lstrip().rstrip(',').split('#', 1)
    return content[0]


def walltimetoseconds(i):
    ctt = i.split('-')
    if len(ctt) < 2:
        pt = datetime.datetime.strptime(ctt[0], '%H:%M:%S')
        total_seconds = pt.second + pt.minute * 60 + pt.hour * 3600
    else:
        total_seconds = eval(ctt[0]) * 24 * 3600
        content = ctt[1].split(':')
        for t in range(3):
            if content[t][0] == '0':
                content[t] = content[t][1:]
        total_seconds += eval(content[2]) + eval(content[1]) * 60 + eval(content[0]) * 3600
    return total_seconds


def distance(a, b):
    t = 0
    for i in range(3):
        t += (b[i] - a[i]) * (b[i] - a[i])
    return np.sqrt(t)


def writeTk(Atoms):
    des = open('T.txt', 'w')
    a = vec(Atoms[4], Atoms[12])
    for i in range(7):
        b = vec(Atoms[i * 9 + 4], Atoms[i * 9 + 12])
        des.write('%d %f %f\n' % (i, k(a, b), distance(Atoms[4], Atoms[i * 9 + 4])))
    des.close()


def k(a, b):
    return np.dot(a, b) / leng(a) / leng(b)


def vec(a, b):
    return [b[i] - a[i] for i in range(3)]


def leng(a):
    t = 0
    for i in range(3):
        t += a[i] * a[i]
    return np.sqrt(t)


def box(ifile, logfile, ofile):
    src = open(ifile, 'r')
    des = open(ofile, 'w')
    Box = readlog(logfile)
    rep = open('replicate.txt', 'w')
    x = int(30 / Box / np.sqrt(3)) + 1
    y = int(30 / Box) + 1
    rep.write('replicate %d %d 1\n' % (x, y))
    rep.close()
    for line in src.readlines():
        ln = line.split()
        if 'xlo' in ln:
            xmin = eval(ln[0])
            xmax = eval(ln[1])
            Xmin = (xmin + xmax) / 2 - Box * np.sqrt(3) / 4
            Xmax = (xmin + xmax) / 2 + Box * np.sqrt(3) / 4
            ln[0] = str(Xmin)
            ln[1] = str(Xmax)
            des.write(' '.join(ln) + '\n')
            continue
        if 'ylo' in ln:
            ymin = eval(ln[0])
            ymax = eval(ln[1])
            Ymin = (ymin + ymax) / 2 - Box / 4
            Ymax = (ymin + ymax) / 2 + Box / 4
            ln[0] = str(Ymin)
            ln[1] = str(Ymax)
            des.write(' '.join(ln) + '\n')
            continue
        des.write(line)
    src.close()
    des.close()


def readlog(logfile):
    src = open(logfile, 'r')
    t = 0
    for line in src.readlines():
        ln = line.split()
        if 'Lx' in ln:
            Lx = ln.index('Lx')
            Ly = ln.index('Ly')
            Den = ln.index('Density')
            t = 1
            continue
        if t:
            lx = eval(ln[Lx])
            ly = eval(ln[Ly])
            den = eval(ln[Den])
            break
    L = np.sqrt(den * lx * ly / 1.2 / np.sqrt(3))
    src.close()
    return L


def main(a, b, c):
    src = open('tablett.txt', 'r')
    if (not (os.path.exists('structure'))):
        os.mkdir('structure')
    os.chdir('structure')
    for line in src.readlines():
        ln = line.split()
        if len(ln) == 1:
            tac = ln[0]
            continue
        if (not (os.path.exists(tac))):
            os.mkdir(tac)
        os.chdir(tac)
        NoM = eval(ln[3])
        tit = ln[0] + ln[1]
        if (not (os.path.exists(tit))):
            os.mkdir(tit)
        else:
            os.chdir('..')
            continue
        os.chdir(tit)
        for i in range(10):
            os.mkdir(str(i))
            os.chdir('../../../polymod/pan')
            sentence = 'stereo s1 '
            replace = ('%s 2 ' % ('weight' if tac == 'ata' else 'pattern')) + ln[0] + (
                ' 0.5 ' if tac == 'ata' else ' ') + ln[1] + (' 0.5' if tac == 'ata' else '')
            copy_file('incha.txt', 'tmpcha.txt', sentence, replace)
            sentence = 'monomers '
            copy_file('tmpcha.txt', 'cha.txt', sentence, str(NoM + 2))
            os.system('../bin/latch cha.txt')
            os.system('python main_test.py')
            while not (validate('graphene.data') and valibond('graphene.data')):
                os.system('../bin/latch cha.txt')
                os.system('python main_test.py')
            os.system('cp graphene.data ../../structure/%s/%s/%d' % (tac, tit, i))
            os.chdir('../../structure/%s/%s/%d' % (tac, tit, i))
            os.system('cp ../../../../input/change.in .')
            os.system('cp ../../../../input/kaka.m .')
            change_kaka('graphene.data', 'polymer_relax.data')
            os.system('cat change.in discommend.txt > chichi.in')
            os.system('lmp < chichi.in')
            change_box('new.data', 'polymer_relax.data')
            os.system('cat ../../../../input/title.in group.txt ../../../../input/link.in > link.in')
            os.system('lmp < link.in')
            os.system('cp ../../../../input/chirality.in .')
            os.system('cp mini0.data polymer_relax.data')
            os.system('lmp < chirality.in')
            box('mini0.data', 'log.lammps', 'polymer_relax.data')
            os.system('cp ../../../../input/glance1.txt .')
            os.system('cat chirality.in replicate.txt glance1.txt > new.in')
            os.system('lmp < new.in')
            os.system("mkdir Qeq")
            for j in ['Qeq']:
                os.system('cp ../../../../input/cha.txt %s' % (j))
                os.system('cp ../../../../input/%snvt.in %s' % (j, j))
                os.system('cp ../../../../input/%snpt.in %s' % (j, j))
                os.system('cp ../../../../input/nptsub.py %s' % (j))
                os.chdir(j)
                os.system('cat ../../../../../input/title.in cha.txt %snvt.in > nvt.in' % j)
                os.system('cp ../glance1.data polymer_relax.data')
                if j == 'Qeq':
                    os.system('cp ../../../../../input/param.qeq.pan .')
                os.system('cat ../../../../../input/title.in cha.txt %snpt.in > npt.in' % j)
                #write_sublmp(tac,tit,i,j,'nvt.in')
                write_sublmp(7200, "%s%d%s%s" % ('nvt', i, tit, tac), 'nvt.in', tmp=0, te='')
                os.system('sbatch sub.sh')
                os.chdir('..')
            os.chdir('..')
        os.chdir('../..')
    src.close()
