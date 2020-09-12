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

from .types import AtomType
from .utils import DEG2RAD
from .config import REAL_MAX

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
    # getPosition()
    # ----------------------------------------------------------------------------
    # Result: store the position of the index-th atom in *pos; calls choke() if
    # index is out of range
    # ============================================================================
    def getPosition(self, index, pos):
        if index < 0 or index >= self.num_entries:
            raise ValueError("getPosition(): Position index %d is out of range (%d)" % (index, self.num_entries))
        if index < self.num_positions and self.positions[index][0] > NO_POSITION:
            pos = self.positions[index]
        else:
            a = np.zeros(3)
            b = np.zeros(3)
            c = np.zeros(3)
            n1 = np.zeros(3)
            n2 = np.zeros(3)
            bond_index_pos = np.zeros(3)
            angle_atom_pos = np.zeros(3)
            torsion_atom_pos = np.zeros(3)
            ze = self.entries[index]

            self.getPosition(ze.bond_index, bond_index_pos)
            bond_length = ze.bond_length
            if -1 == ze.angle_index:
                # Atom 1
                pos[0] = bond_length
                pos[1] = 0.0
                pos[2] = 0.0
            else:
                self.getPosition(ze.angle_index, angle_atom_pos)
                bond_angle = DEG2RAD(ze.bond_angle)
                if -1 == ze.dihedral_index:
                    # Atom 2
                    torsion_atom_pos[0] = 0.0
                    torsion_atom_pos[1] = 1.0
                    torsion_atom_pos[2] = 0.0
                    torsion_angle = DEG2RAD(90.0)
                else:
                    # All other atoms
                    self.getPosition(ze.dihedral_index, torsion_atom_pos)
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
    for i in range(zm.num_positions):
        zm.positions[i] = np.zeros(3)
        zm.positions[i][0] = NO_POSITION
        zm.positions[i][1] = NO_POSITION
        zm.positions[i][2] = NO_POSITION
    return zm
