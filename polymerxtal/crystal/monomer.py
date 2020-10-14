"""
Functions for calculating monomer properties.
"""

from polymerxtal.io import run_polymod, read_latchout
from polymerxtal.polymod import readPDB

from .molecule import build_bond_list


# A bond between two monomer atoms, not represented by the z-matrix
class Backbone:
    def __init__(self, name, path, head_index, tail_index, tail_torsion):
        self.name = name
        self.path = path
        self.head_index = head_index
        self.tail_index = tail_index
        self.tail_torsion = tail_torsion


def get_tail_torsion(mpath, head, tail):
    logfile = 'latch.log'
    infile = 'latch.in'
    des = open(infile, 'w')
    des.write('bond_scale  1.2\n')
    des.write('monomer m1 %s %d %d\n' % (mpath, head, tail))
    des.write('torsion all fixed \n')
    des.write('chains 1\n')
    des.write('monomers 1\n')
    des.write('system_min 0 0 0\n')
    des.write('system_max 100 100 100\n')
    des.write('grid_size 5\n')
    des.write('stereo s pattern 1 m1\n')
    des.write('chain_stereo 1 s 1\n')
    des.write('sample none\n')
    des.write('energy none\n')
    des.write('write intermediate\n')
    des.write('log_file %s\n' % logfile)
    des.close()
    run_polymod(infile)
    latch_info = read_latchout(logfile)
    tail_torsion = eval(latch_info['monomer'][1]['info'][3][-1])
    return tail_torsion


def identify_backbone_path(monomer_path, center_atoms, side_atom=0):
    h = readPDB(monomer_path)
    bonds = build_bond_list(h.pos)

    center_head = center_atoms[0]
    center_tail = center_atoms[-1]
    if side_atom and ((center_tail - 1, side_atom - 1) not in bonds) and ((side_atom - 1, center_tail - 1)
                                                                          not in bonds):
        raise ValueError('side atom must connect to tail backbone atom')

    head_atoms = []
    tail_atoms = []

    for key in bonds:
        atom1, atom2 = key
        if side_atom and (side_atom - 1 in [atom1, atom2]):
            continue
        elif center_head - 1 == atom1:
            if not center_tail - 1 == atom2:
                head_atoms.append(atom2 + 1)
        elif center_head - 1 == atom2:
            if not center_tail - 1 == atom1:
                head_atoms.append(atom1 + 1)
        elif center_tail - 1 == atom1:
            if not center_head - 1 == atom2:
                tail_atoms.append(atom2 + 1)
        elif center_tail - 1 == atom2:
            if not center_head - 1 == atom1:
                tail_atoms.append(atom1 + 1)

    backbones = {}
    tail_flag = 0

    for tail in tail_atoms:
        tail_flag += 1
        for head in head_atoms:
            tail_torsion = get_tail_torsion(monomer_path, head, tail)
            while tail_torsion < 0:
                tail_torsion += 360
            while tail_torsion >= 360:
                tail_torsion -= 360
            tail_torsion = 60 if 0 <= tail_torsion < 120 else 180 if 120 <= tail_torsion < 240 else 300
            m_name = 'm' + str(tail_flag) + 'T' + str(tail_torsion)
            backbones[m_name] = Backbone(m_name, monomer_path, head, tail, tail_torsion)

    return backbones
