# ============================================================================
# chain.py -- Chain functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .zmatrix import ZMatrix, ZEntry, MIN_POSITIONS
from .stereo import Stereo
from .monomer import Bond, Monomer
from .stdio import FILE
from .vector import Vector
from .utils import foldPosition
from .element import getElementName
from .params import Params


class Chain:
    def __init__(self):
        self.zm = ZMatrix()
        self.stereo = Stereo()
        self.num_monomers = 0
        self.curr_monomer = 0
        self.curr_atom = 0
        self.tail_index = -1
        self.init_domain = 0
        self.domain = 0
        self.dead = 0
        self.i_monomer = {}  # indices of monomers within the Chain
        self.torsion_count = {}
        self.length = {}  # end to end, as each monomer is added
        self.weight = {}  # as each monomer is added
        self.mass = 0.
        self.extra_bonds = Bond()
        self.extra_bonds.create()

    # ============================================================================
    # destroyChain()
    # ----------------------------------------------------------------------------
    # Result: free a Chain
    # ============================================================================

    def __del__(self):
        if self:
            if self.zm:
                del self.zm
                self.zm = ZMatrix()
            self.i_monomer = {}
            self.length = {}
            self.weight = {}
            del self.extra_bonds
            self.extra_bonds = Bond()
            self.torsion_count = {}

    # ============================================================================
    # addMonomer()
    # ----------------------------------------------------------------------------
    # Result: add a new Monomer of type m to the Chain c
    # ============================================================================

    def addMonomer(self, m, bond_length, update_bonds):
        na = self.curr_atom
        czm = self.zm
        mzm = m.zm
        pos = Vector()

        self.curr_monomer += 1
        self.i_monomer[c.curr_monomer] = na
        if 0 == na:
            # First monomer in c
            for i in range(mzm.num_entries):
                czm.entries[i] = mzm.entries[i]  # struct copy
            self.curr_atom = mzm.num_entrie
            offset = 0
            self.mass += m.head_mass
        else:
            old_tail_index = self.tail_index
            cze = czm.entries[old_tail_index]
            mze = ZEntry()

            offset = na - 2
            # Replace current tail of c with new monomer head neighbor
            cze.type = mzm.entries[1].type  # struct copy
            cze.bond_length = bond_length
            if self.zm.num_positions > MIN_POSITIONS:
                # Storing positions -- update previous tail atom position
                czm.clearPosition(old_tail_index)
                czm.getPosition(old_tail_index, pos)
                czm.setPosition(old_tail_index, pos)
            # Add remaining atoms in monomer
            for i in range(2, mzm.num_entries):
                self.curr_atom += 1
                cze = czm.entries[c.curr_atom]
                mze = mzm.entries[i]
                cze = mze  # struct copy
                cze.bond_index = shiftIndex(czm, mze.bond_index, old_tail_index, offset)
                cze.angle_index = shiftIndex(czm, mze.angle_index, old_tail_index, offset)
                cze.dihedral_index = shiftIndex(czm, mze.dihedral_index, old_tail_index, offset)
        self.tail_index = offset + m.num_bb - 1
        self.mass += m.central_mass
        if self.curr_monomer == self.num_monomers:  # last monomer
            self.mass += m.tail_mass
        if update_bonds:
            blist = m.extra_bonds
            b = Bond()
            b.create()

            while blist and hasattr(blist, 'next'):
                b = createBond(blist.index1 + offset, blist.index2 + offset)
                b.next = self.extra_bonds
                self.extra_bonds = b
                blist = blist.next

    # ============================================================================
    # writeChainPDB()
    # ----------------------------------------------------------------------------
    # Result: write a PDB representation of the Chain c to the FILE f using
    # unwrapped coordinates if fold_pos is 0 or wrapped coordinates if fold_pos
    # is 1; increment *i_atom as each atom is written, and roll *i_atom back to 1
    # after 99999
    # ============================================================================

    def writeChainPDB(self, i_atom, fold_pos, f, p):
        pos = Vector()

        RESIDUE_NAME = "UNK"
        CHAIN_ID = "A"
        RESIDUE_SEQ = 1

        for i in range(self.curr_atom):
            self.zm.getPosition(i, pos)
            if fold_pos:
                foldPosition(pos, p.system_min, p.system_max, p.system_size)
            #        1     7   13  18     23     31   39   47
            f.write("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f\n" %
                    (i_atom, getElementName(self.zm.entries[i].type.element_index), RESIDUE_NAME, CHAIN_ID,
                     RESIDUE_SEQ, pos.x, pos.y, pos.z))
            i_atom += 1
            if i_atom > 99999:
                i_atom = 1
        i_atom += 1
        #  1     7
        f.write("TER   %5d\n" % i_atom)


# ============================================================================
# createChain()
# ----------------------------------------------------------------------------
# Result: return a newly allocated, initialized Chain that will contain
# num_monomers repeated units and num_atoms total atoms
# ============================================================================


def createChain(stereo, num_monomers, num_atoms, store_positions):
    c = Chain()

    c.zm = createZMatrix(num_atoms, store_positions)
    c.stereo = stereo
    for i in range(num_monomers):
        c.i_monomer[i] = 0
        c.length[i] = 0.
        c.weight[i] = 0.
    c.num_monomers = num_monomers
    for i in range(360):
        c.torsion_count[i] = 0.
    c.extra_bonds = Bond()
    c.extra_bonds.create()
    c.init_domain = -1
    return c


# ============================================================================
# shiftIndex()
# ----------------------------------------------------------------------------
# Result: return the atom index of an atom in a ZMatrix added to an existing
# ZMatrix
# ============================================================================


def shiftIndex(zm, index, old_tail, offset):
    retval = -1

    if index == -1:
        retval = zm.entries[zm.entries[old_tail].bond_index].bond_index
    elif index == 0:
        retval = zm.entries[old_tail].bond_index
    elif index == 1:
        retval = old_tail
    else:
        retval = index + offset
    return retval


# ============================================================================
# writeInternalRotationPDB()
# ----------------------------------------------------------------------------
# Result: write a series of PDB files representing rotation of a backbone atom
# in a monomer of type m about the bond joining it to its neighbor
# ============================================================================


def writeInternalRotationPDB(m, index, step, basename):
    c = createChain(Stereo(), 1, m.num_atoms, 0)
    leng = len(basename) + 8  # basename_xxx.pdb
    f = FILE()

    c.addMonomer(m, 0.0, 0)
    for i in range(0, 360, step):
        c.zm.entries[index].torsion_angle = float(i)
        path = "%s_%3d.pdb" % (basename, i)
        f = open(path, "w")
        n = 1
        c.writeChainPDB(n, 0, f, Params())
        f.close()
    del c
    path = ''
