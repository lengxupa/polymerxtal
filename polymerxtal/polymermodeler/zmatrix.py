# ============================================================================
# zmatrix.py -- ZEntry, ZMatrix structs; prototypes for functions used in
# latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

from .config import REAL_MAX, MAX_BONDS
from .element import getElementName
from .stdio import FILE
from .types import AtomType
from .utils import DEG2RAD

NO_POSITION = -REAL_MAX

# Minimum number of stored positions
MIN_POSITIONS = 3


# Internal coordinate representation of a single atom
class ZEntry:
    def __init__(self):
        self.type = AtomType()
        self.bond_index = -1
        self.angle_index = -1
        self.dihedral_index = -1
        self.bond_length = 0.  # Angstroms
        self.bond_angle = 0.  # degrees
        self.torsion_angle = 0.  # degrees


class ZMatrix:
    def __init__(self):
        self.num_entries = 0  # size of entries[]
        self.num_positions = 0
        self.entries = {}
        self.positions = {}

    # ============================================================================
    # setPosition()
    # ----------------------------------------------------------------------------
    # Result: set the position of the index-th atom; calls choke() is index is out
    # of range
    # ============================================================================
    def setPosition(self, index, pos):
        if index < 0 or index >= self.num_positions:
            raise ValueError("setPosition(): Position index %d is out of range" % index)
        self.positions[index] = pos

    # ============================================================================
    # clearPosition()
    # ----------------------------------------------------------------------------
    # Result: set the position of the index-th atom to NO_POSITION; calls choke()
    # if index is out of range
    # ============================================================================
    def clearPosition(self, index):
        if index < 0 or index >= self.num_positions:
            raise ValueError("clearPosition(): Position index %d is out of range" % index)
        self.positions[index][0] = NO_POSITION
        # Only x is tested in getPosition()
        # zm->positions[index].y = NO_POSITION;
        # zm->positions[index].z = NO_POSITION;

    # ============================================================================
    # isBonded()
    # ----------------------------------------------------------------------------
    # Result: return 1 if entries n and m are separated by at most max_bonds
    # bonds in zm, else return 0
    # ============================================================================
    def isBonded(self, n, m, max_bonds):
        m_bonds = {}
        n_bonds = {}
        for i in range(MAX_BONDS):
            m_bonds[i] = 0
            n_bonds[i] = 0

        m_bonds[0] = self.entries[m].bond_index
        i = 1
        while i < max_bonds and m_bonds[i - 1] > -1:
            m_bonds[i] = self.entries[m_bonds[i - 1]].bond_index
            if n == m_bonds[i]:
                return 1
            i += 1
        m_count = i

        n_bonds[0] = self.entries[n].bond_index
        i = 1
        while i < max_bonds and n_bonds[i - 1] > -1:
            n_bonds[i] = self.entries[n_bonds[i - 1]].bond_index
            if m == n_bonds[i]:
                return 1
            i += 1
        n_count = i

        for i in range(m_count):
            common = -1
            j = 0
            while -1 == common and j < n_count:
                if m_bonds[i] == n_bonds[j]:
                    common = i
                else:
                    j += 1
            if common > -1 and (i + j + 2) <= max_bonds:
                return 1
        return 0

    # ============================================================================
    # getPosition()
    # ----------------------------------------------------------------------------
    # Result: store the position of the index-th atom in *pos; calls choke() if
    # index is out of range
    # ============================================================================
    def getPosition(self, index):
        if index < 0 or index >= self.num_entries:
            raise ValueError("getPosition(): Position index %d is out of range (%d)" % (index, self.num_entries))
        if index < self.num_positions and self.positions[index][0] > NO_POSITION:
            pos = self.positions[index]
        else:
            torsion_atom_pos = np.zeros(3)
            ze = self.entries[index]

            bond_index_pos=self.getPosition(ze.bond_index)
            bond_length = ze.bond_length
            if -1 == ze.angle_index:
                # Atom 1
                pos[0] = bond_length
                pos[1] = 0.0
                pos[2] = 0.0
            else:
                angle_atom_pos=self.getPosition(ze.angle_index)
                bond_angle = DEG2RAD(ze.bond_angle)
                if -1 == ze.dihedral_index:
                    # Atom 2
                    torsion_atom_pos[0] = 0.0
                    torsion_atom_pos[1] = 1.0
                    torsion_atom_pos[2] = 0.0
                    torsion_angle = DEG2RAD(90.0)
                else:
                    # All other atoms
                    torsion_atom_pos=self.getPosition(ze.dihedral_index)
                    torsion_angle = DEG2RAD(ze.torsion_angle)
                a = bond_index_pos - angle_atom_pos
                b = bond_index_pos - torsion_atom_pos
                n1 = np.cross(a, b)
                n1 /= np.linalg.norm(n1)
                n2 = np.cross(a, n1)
                n2 /= np.linalg.norm(n2)
                n1 *= -np.sin(torsion_angle)
                n2 *= np.cos(torsion_angle)
                c = n1 + n2
                c /= np.linalg.norm(c)
                c *= bond_length * np.sin(bond_angle)
                b = bond_index_pos + c
                a /= np.linalg.norm(a)
                a *= bond_length * np.cos(bond_angle)
                pos = b - a
        return pos

    # ============================================================================
    # writeZMatrix()
    # ----------------------------------------------------------------------------
    # Result: write a ZMatrix to a FILE
    # ============================================================================
    def writeZMatrix(self, f, Dreiding_type):
        for i in range(self.num_entries):
            if Dreiding_type:
                self.entries[i].type.writeAtomTypeDreiding(f)
            else:
                f.printf("%2.2s" % getElementName(self.entries[i].type.element_index))
            if 0 == i:
                f.printf("\n")
                continue
            f.printf("  %2d  %5.2f" % (self.entries[i].bond_index + 1, self.entries[i].bond_length))
            if 1 == i:
                f.printf("\n")
                continue
            f.printf("  %2d  %6.2f" % (self.entries[i].angle_index + 1, self.entries[i].bond_angle))
            if 2 == i:
                f.printf("\n")
                continue
            f.printf("  %2d  %7.2f" % (self.entries[i].dihedral_index + 1, self.entries[i].torsion_angle))


# ============================================================================
# createZMatrix()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized ZMatrix; calls
# choke() if allocation fails
# ============================================================================
def createZMatrix(num_entries, store_positions):
    zm = ZMatrix()

    zm.num_entries = num_entries
    for i in range(num_entries):
        zm.entries[i] = ZEntry()
        zm.entries[i].type.element_index = -1
        zm.entries[i].type.num_bonds = 0
        zm.entries[i].bond_index = -1
        zm.entries[i].angle_index = -1
        zm.entries[i].dihedral_index = -1
        zm.entries[i].bond_length = 0.0
        zm.entries[i].bond_angle = 0.0
        zm.entries[i].torsion_angle = 0.0
    zm.num_positions = num_entries if store_positions else MIN_POSITIONS
    for i in range(zm.num_positions):
        zm.positions[i] = np.zeros(3)
        zm.positions[i][0] = NO_POSITION
        zm.positions[i][1] = NO_POSITION
        zm.positions[i][2] = NO_POSITION
    return zm
