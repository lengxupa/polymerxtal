# ============================================================================
# monomer.py -- Monomer typedef, prototypes for functions used in latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from enum import Enum
import numpy as np

from polymerxtal import calculate_angle

from .config import REAL_MAX
from .element import getElementIndex, getElementName, getElementMass, getElementR0
from .scan import *
from .stdio import FILE, fgets
from .types import *
from .utils import RAD2DEG, kB, openFile
from .zmatrix import ZMatrix, createZMatrix, ZEntry


# A bond between two monomer atoms, not represented by the z-matrix
class Bond:
    def __init__(self):
        self.index1 = -1
        self.index2 = -1
        self.visited = 0

    def create(self):
        self.next = Bond()


# Restrictions on allowed torsion (dihedral) angles for backbone atoms
class Torsion(Enum):
    TORSION_FIXED = 0  # Torsion may not change
    TORSION_FREE = 1  # Torsion may assume any value with equal probability
    TORSION_ENERGY = 2  # Torsion probability determined by E(phi)
    TORSION_ENERGY_CALC = 3  # E(phi) calculated internally


# Temporary structure used to hold atomic positions and atom types while
# creating and filling a Monomer structure
class Holder:
    def __init__(self):
        self.pos = {}
        self.el_names = {}
        self.atom_types = {}
        self.num_atoms = 0


# Used only in this file; store name-value pairs when parsing zmatrix file
class ZVar:
    def __init__(self):
        self.name = ''
        self.value = ''

    def create(self):
        self.next = ZVar()


# File scope
z_variables = ZVar()


# ============================================================================
# createHolder()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Holder; calls
# choke() on allocation failure
# ============================================================================
def createHolder(num_atoms):
    h = Holder()

    for i in range(num_atoms):
        h.pos[i] = np.zeros(3)
        h.el_names[i] = ''
        h.atom_types[i] = AtomType()
    h.num_atoms = num_atoms
    return h


# ============================================================================
# addZVar()
# ----------------------------------------------------------------------------
# Result: allocates a new ZVar, along with name and value fields, and adds it
# to z_variables list
# ============================================================================
def addZVar(name, value):
    zv = ZVar()
    zv.create()

    zv.name = name
    zv.value = value

    #printf("variable: %s = %s\n", zv->name, zv->value);  fflush(stdout);

    zv.next = z_variables
    z_variables = zv


# ============================================================================
# getRealValue()
# ----------------------------------------------------------------------------
# Result: return the Real value associated with str -- if str is the name of a
# ZVar in the z_variables list, the value of that variable converted to Real
# (with possible recursion, if the value is also the name of another ZVar), or
# str converted to Real
# ============================================================================
def getRealValue(string):
    zv = z_variables

    while zv and hasattr(av, 'next'):
        if zv.name == string:
            return getRealValue(zv.value)
        zv = zv.next
    return float(string)


# ============================================================================
# cleanupZVars()
# ----------------------------------------------------------------------------
# Result: frees all ZVars in the z_variables list; sets z_variables to NULL
# ============================================================================
def cleanupZVars():

    zv1 = z_variables
    while zv1 and hasattr(av1, 'next'):
        zv2 = zv1.next
        zv1.next = ZVars()
        zv1.name = ''
        zv1.value = ''
        del zv1
        zv1 = zv2
    z_variables = ZVars()


# ============================================================================
# readXYZ()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly created Holder that contains the element
# names and atomic positions read from input file
# ============================================================================
def readXYZ(s):
    h = Holder()

    # Expected format of the XYZ file:
    # N                  # line 1 - number of atoms (N)
    # monomer name       # line 2 - arbitrary string (ignored)
    # element  x  y  z   # line 3 - atom 1
    # ...
    # element  x  y  z   # line N+2 - atom N
    #
    # element: string e.g. "C", "H", "O"
    # coordinates: x y z in Angstroms

    n = s.getIntToken()
    h = createHolder(n)
    s.getToken()
    while 3 > s.lineno:
        s.getToken()
    # s->tokstr points to the first atom type, on line 3
    for i in range(h.num_atoms):
        h.el_names[i] = s.tokstr
        h.pos[i][0] = s.getRealToken()
        h.pos[i][1] = s.getRealToken()
        h.pos[i][2] = s.getRealToken()
        s.getToken()
    return h


# ============================================================================
# readZM()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly created Holder that contains the element
# names and positions (with the first atom at the origin); creates and fills
# *zmptr with the input z-matrix, except for the atom types; creates and
# destroys ZVars
# ============================================================================
def readZM(s, zmptr):
    h = {}
    TOKLEN = 256
    v0 = np.zeros(3)

    # Expected format of the zmatrix file:
    # N  # not used
    # type
    # type bond_index bond_length
    # type bond_index bond_length angle_index bond_angle
    # type bond_index bond_length angle_index bond_angle dih_index torsion_angle
    # ...
    # var = value
    # ...
    #
    # type:   element string e.g. "C", "H", "O"
    # length: Angstroms
    # angle:  degrees
    #
    # Length and angle entries can be a variable name.
    #
    # This method reads the file twice; see comments below in findBonds()
    # for why this is not a problem.

    if TOK_INT == s.getToken():
        i = 0
        skip_N = 1
    else:
        i = 1
        skip_N = 0
    # First pass: count atoms, store variables
    # Count tokens until the first '='
    prevtok = "%s" % s.tokstr
    tokval = s.getToken()
    i += 1
    while (tokval != TOK_EOF) and (tokval != TOK_EQUAL):
        # Count tokens until the first '='
        prevtok = "%s" % s.tokstr
        tokval = s.getToken()
        i += 1
    i -= 1  # ignore current token: EOF or '='
    if TOK_EQUAL == tokval:
        i -= 1  # don't count name of first variable
    n = 0
    if i > 0:
        i -= 1  # First atom
        n += 1
    if i > 0:
        i -= 3  # Second atom
        n += 1
    if i > 0:
        i -= 5  # Third atom
        n += 1
    while i > 0:
        # Atoms 4 through N
        i -= 7
        n += 1
    if i != 0:
        raise TypeError("Invalid zmatrix: missing information")
    while TOK_EQUAL == tokval:
        # "name = value" pairs follow; variable name is in prevtok
        s.getToken()
        addZVar(prevtok, s.tokstr)
        s.getToken()
        prevtok = "%s" % s.tokstr
        tokval = s.getToken()

    # Second pass: fill zm
    s.resetScanner()
    if skip_N:
        s.getToken()
    h = createHolder(n)
    zmptr = createZMatrix(n, 0)  # don't store positions
    i = 0
    while i < n:
        s.getToken()
        h.el_names[i] = s.tokstr
        if 0 == i:
            # End first atom
            i += 1
            continue
        # Bond index, length
        zmptr.entries[i].bond_index = s.getIntToken() - 1
        s.getToken()
        zmptr.entries[i].bond_length = getRealValue(s.tokstr)
        if 1 == i:
            # End second atom
            i += 1
            continue
        # Angle index, bond angle
        zmptr.entries[i].angle_index = s.getIntToken() - 1
        s.getToken()
        zmptr.entries[i].bond_angle = getRealValue(s.tokstr)
        if 2 == i:
            # End third atom
            i += 1
            continue
        # Dihedral index, torsion angle
        zmptr.entries[i].dihedral_index = s.getIntToken() - 1
        s.getToken()
        zmptr.entries[i].torsion_angle = getRealValue(s.tokstr)
        i += 1

    # Now, fill h->pos[] from zm internal coordinates
    zmptr.setPosition(0, v0)
    for i in range(h.num_atoms):
        h.pos[i] = zmptr.getPosition(i, h.pos[i])

    #for (i = 0; i < h->num_atoms; i++)
    #{
    #   printf("Atom %d: (%g, %g, %g)\n", i+1, h->pos[i].x, h->pos[i].y,
    #       h->pos[i].z);
    #}
    #fflush(stdout);

    # Cleanup
    cleanupZVars()
    return h


# ============================================================================
# readPDB()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly created Holder that contains the element
# names and atomic positions read from input file; calls choke() if the file
# can't be opened
# ============================================================================
def readPDB(path):
    LINELEN = 81

    # PDB file format: only use the ATOM/HETATM lines
    # Columns        Type            Information
    #  1 -  6        Record name     "ATOM  " or "HETATM"
    #  7 - 11        Integer         Atom serial number
    # 13 - 16        Atom            Atom name
    # 17             Character       Alternate location indicator
    # 18 - 20        Residue name    Residue name
    # 22             Character       Chain identifier
    # 23 - 26        Integer         Residue sequence number
    # 27             AChar           Code for insertion of residues
    # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms
    # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms
    # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms
    # 55 - 60        Real(6.2)       Occupancy
    # 61 - 66        Real(6.2)       Temperature factor (Default = 0.0)
    # 73 - 76        LString(4)      Segment identifier, left-justified
    # 77 - 78        LString(2)      Element symbol, right-justified
    # 79 - 80        LString(2)      Charge on the atom
    #
    # Because many of these fields may not be present, the token parsing
    # methods used in the other two read functions are not as reliable as
    # grabbing one line at a time and examining the specific columns.
    #
    # This method reads the file twice; see comments below in findBonds()
    # for why this is not a problem.

    f = openFile(path, "r")
    n = 0
    line = fgets(LINELEN, f)
    while line:
        if (line[:6] == "ATOM  ") or (line[:6] == "HETATM"):
            n += 1
        line = fgets(LINELEN, f)
    h = createHolder(n)

    f.rewind()
    i = 0
    line = fgets(LINELEN, f)
    while line:
        leng = len(line)
        if (line[:6] == "ATOM  ") or (line[:6] == "HETATM"):
            h.pos[i][0] = float(line[30:38])
            h.pos[i][1] = float(line[38:46])
            h.pos[i][2] = float(line[46:54])
            # Check element symbol first...
            if leng > 77 and (not line[77].isspace()):
                # line[78] = '\0';
                p = 77
                if not line[76].isspace():
                    p -= 1
            else:
                # ...otherwise try to grab atom name, which might have extra
                # information e.g. hybridization

                p = 12
                while line[p].isspace() and p < 16:
                    p += 1
                q = p + 1
                while not line[q].isspace():
                    q += 1
            h.el_names[i] = line[p:q]
            i += 1
        line = fgets(LINELEN, f)
    f.fclose()
    return h


# ============================================================================
# createBond()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Bond; calls
# choke() if allocation fails
# ============================================================================
def createBond(index1, index2):
    b = Bond()

    b.index1 = index1
    b.index2 = index2
    b.visited = 0
    b.next = Bond()
    return b


# ============================================================================
# clearBondVisits()
# ----------------------------------------------------------------------------
# Result: mark every Bond in the list blist un-visited
# ============================================================================
def clearBondVisits(blist):
    b = blist

    while b and hasattr(b, 'next'):
        b.visited = 0
        b = b.next


# ============================================================================
# findBonds()
# ----------------------------------------------------------------------------
# Result: return a list of Bonds between atoms in h; calls choke() if
# h->el_names[] contains an unknown element
# ============================================================================
def findBonds(h, eq_bond_scale):
    blist = Bond()
    bcurr = Bond()

    # This is an O(N^2) operation, but...
    #   1. The number of atoms, N, in a monomer is small
    #   2. This is only done once, during intialization
    for i in range(h.num_atoms):
        el_index = getElementIndex(h.el_names[i])
        dj = getElementR0(el_index)
        for j in range(i + 1, h.num_atoms):
            el_index = getElementIndex(h.el_names[j])
            di = getElementR0(el_index)
            d = np.linalg.norm(h.pos[i] - h.pos[j])
            if d < (di + dj) * eq_bond_scale:
                # Bond
                b = createBond(i, j)
                if blist and hasattr(blist, 'next') and bcurr and hasattr(bcurr, 'next'):
                    bcurr.next = b
                else:
                    blist = b
                bcurr = b

                # printf("Created bond between %d and %d\n", i+1, j+1);

            # else:
            #    printf("No bond between %d and %d\n", i+1, j+1);
            # printf("di = %g, dj = %g, d = %g\n\n", di, dj, d);
    return blist


def FIND_BOND(b, ndx):
    while b and hasattr(b, 'next') and ndx != b.index1 and ndx != b.index2:
        b = b.next
    return b


def OTHER_INDEX(b, ndx):
    return b.index2 if ndx == b.index1 else b.index1


# ============================================================================
# findRing()
# ----------------------------------------------------------------------------
# Result: return 1 if a closed loop, starting and ending at target, can be
# identified, else return 0; recurses
# ============================================================================
def findRing(target, prev, current, blist):
    b = blist

    while b and hasattr(b, 'next'):
        b = FIND_BOND(b, current)
        if b and hasattr(b, 'next'):
            n = OTHER_INDEX(b, current)
            if n != prev and (not b.visited):
                b.visited += 1
                if target == n:
                    return 1  # bingo -- a closed loop
                else:
                    if findRing(target, current, n, blist):
                        return 1
            b = b.next
    return 0


# ============================================================================
# setAtomTypes()
# ----------------------------------------------------------------------------
# Result: fills h->atom_types and registers AtomTypes; calls choke() if an
# unknown element name is found
# ============================================================================
def setAtomTypes(h, blist):
    H_index = getElementIndex("H")

    for i in range(h.num_atoms):
        at = AtomType()
        # Count the number of bonds which involve atom i
        at.element_index = getElementIndex(h.el_names[i])
        at.num_bonds = 0
        ring = 0
        b = blist
        while b and hasattr(b, 'next'):
            if i == b.index1 or i == b.index2:
                at.num_bonds += 1
                clearBondVisits(blist)
                b.visited += 1
                if (not ring) and findRing(i, i, OTHER_INDEX(b, i), blist):
                    ring = 1
            b = b.next
        if ring:  # X_R
            # printf("atom %d is in a resonant bond\n", i+1); fflush(stdout);
            at.num_bonds *= -1
        h.atom_types[i] = at
    # Identify hydrogen bonds and register new types: do this as a separate
    # loop so we know what atom type is bonded to each H

    for i in range(h.num_atoms):

        if H_index == h.atom_types[i].element_index:
            b = blist
            b = FIND_BOND(b, i)  # only 1 bond
            name = getElementName(h.atom_types[OTHER_INDEX(b, i)].element_index)
            if "N" == name or "O" == name or "F" == name:
                # printf("atom %d is in a H-bond\n", i+1); fflush(stdout);
                h.atom_types[i].num_bonds = -1

        if -1 == getAtomTypeIndex(h.atom_types[i]):  # new type
            addAtomType(h.atom_types[i])


# ============================================================================
# setBondTypes()
# ----------------------------------------------------------------------------
# Result: registers known bond types
# ============================================================================
def setBondTypes(h, blist):
    at1 = AtomType()
    at2 = AtomType()

    b = blist
    while b and hasattr(b, 'next'):
        at1 = h.atom_types[b.index1]
        at2 = h.atom_types[b.index2]
        if -1 == getBondTypeIndex(at1, at2):
            addBondType(at1, at2)
        b = b.next


TARGET_SQ_DIST = 0.01

ANGLE_DIFF = 60.0


# Entry in a linked list
class Monomer:
    def __init__(self):
        self.name = ''
        self.num_atoms = 0
        self.num_extra_bonds = 0
        self.num_bb = 0  # Number of backbone atoms; head: 0;  tail: num_bb-1
        self.selection_count = 0
        self.zm = {}
        self.extra_bonds = Bond()
        self.extra_bonds.create()
        self.torsions = {}
        self.torsion_angles = {}
        self.torsion_energies = {}
        self.torsion_probs = {}
        self.torsion_prob_min = {}
        self.head_mass = 0.
        self.tail_mass = 0.
        self.central_mass = 0.  # all atoms except head and tail
        self.num_torsions = {}

    def create(self):
        self.next = Monomer()

    # ============================================================================
    # fillZMatrix()
    # ----------------------------------------------------------------------------
    # Result: fills m->zm and sets m->num_bb; calls choke() on fatal errors --
    # inability to identify backbone bonds from head to tail, NaN angles,
    # allocation errors
    # ============================================================================
    def fillZMatrix(self, h, blist, head_index, tail_index):
        indices = {}
        rev_indices = {}
        not_bb = {}
        i_not_bb = 0
        v0 = np.zeros(3)

        # Use indices[] to order the atoms, starting with the head, and identify
        # backbone atoms to the tail.  The chain building algorithms depend on
        # having the N backbone atoms in the first N positions in the ZMatrix.
        indices[0] = head_index
        rev_indices[head_index] = 0
        have_backbone = 0
        i = 1
        while not have_backbone:
            b = blist
            d2min = REAL_MAX
            bnear = Bond()
            while b and hasattr(b, 'next'):
                # Visit all bonds involving atom indices[i-1]; find the bonded atom
                # nearest to the tail atom.
                b = FIND_BOND(b, indices[i - 1])
                if b and hasattr(b, 'next'):
                    j = 0
                    while j < i_not_bb:
                        if OTHER_INDEX(b, indices[i - 1]) == not_bb[j]:
                            break
                        j += 1
                    if j == i_not_bb:
                        # atom to which indices[i-1] is bonded has not already been
                        # rejected from the backbone
                        d2 = np.linalg.norm(h.pos[OTHER_INDEX(b, indices[i - 1])] - h.pos[tail_index])**2
                        if d2 < d2min:
                            d2min = d2
                            bnear = b
                    b = b.next
            if bnear and hasattr(bnear, 'next'):
                n = OTHER_INDEX(bnear, indices[i - 1])
                bond_used = 0
                for j in range(i):
                    if indices[j] == n:
                        bond_used += 1
                if not bond_used:
                    indices[i] = n
                    rev_indices[indices[i]] = i
                    self.zm.entries[i].bond_index = i - 1

                    #printf("backbone atom %d: %d\n", i+1, indices[i]+1); fflush(stdout);

                    i += 1
                    if n == tail_index:
                        have_backbone += 1
                        self.num_bb = i

                        #printf("found tail!\n"); fflush(stdout);

                else:
                    # Chosen bond (nearest to tail) is already in the backbone!
                    # Current backbone atom is not truly on the backbone.
                    not_bb[i_not_bb] = indices[i - 1]
                    i_not_bb += 1
                    i -= 1
            else:
                raise TypeError("Unable to fill zmatrix")
        # printf("num_bb = %d\n", m->num_bb); fflush(stdout);

        # Now i is the index of the first non-backbone atom.  Add the remaining
        # atom indices, ensuring that we add indices of atoms that are bonded to
        # indices that are already in the z-matrix
        while i < self.num_atoms:
            for j in range(i):
                # Find a bond involving atom indices[j]
                b = blist
                b = FIND_BOND(b, indices[j])
                while b and hasattr(b, 'next'):
                    k = OTHER_INDEX(b, indices[j])
                    n = 0
                    while n < i:
                        if k == indices[n]:
                            break
                        n += 1
                    if n == i:  # k is not yet in the z-matrix
                        indices[i] = k

                        # printf("indices[%d] = %d\n", i, indices[i]); fflush(stdout);

                        rev_indices[indices[i]] = i
                        self.zm.entries[i].bond_index = j
                        i += 1
                    else:
                        b = b.next
                        b = FIND_BOND(b, indices[j])

        # printf("indices[]: ");
        # for (i = 0; i < m->num_atoms; i++)
        #    printf("%d ", indices[i]);
        # printf("\n");
        # printf("rev_indices[]: ");
        # for (i = 0; i < m->num_atoms; i++)
        #    printf("%d ", rev_indices[i]);
        # printf("\n");
        # printf("Bond indices: ");
        # for (i = 0; i < m->num_atoms; i++)
        #    printf("%d ", m->zm->entries[i].bond_index);
        # printf("\n");
        # fflush(stdout);

        # Now indices[] holds the order of atoms in the ZMatrix, with the backbone
        # atoms first; indices[i] maps the position i in the ZMatrix to the
        # positions in h->pos[].  We also have a reverse mapping, rev_indices[],
        # which maps the atom index j in h->pos[j] to the ZMatrix position, i.

        # Types
        for i in range(self.num_atoms):
            self.zm.entries[i].type = h.atom_types[indices[i]]

        # Bond lengths
        for i in range(1, self.num_atoms):
            j = indices[i]
            k = indices[self.zm.entries[i].bond_index]
            self.zm.entries[i].bond_length = np.linalg.norm(h.pos[j] - h.pos[k])

        # Angle indices, bond angles
        for i in range(2, self.num_atoms):
            j = self.zm.entries[i].bond_index
            k = self.zm.entries[j].bond_index
            self.zm.entries[i].angle_index = k
            v1 = h.pos[indices[i]] - h.pos[indices[j]]
            v2 = h.pos[indices[k]] - h.pos[indices[j]]
            self.zm.entries[i].bond_angle = RAD2DEG(calculate_angle(-v1, np.zeros(3), -v2))
            if np.isnan(self.zm.entries[i].bond_angle):
                raise ValueError("Bond angle for atom %d is NaN" % (i + 1))

        # Dihedral index, torsion angles
        # setPosition(m->zm, 0, &v0);
        self.zm.setPosition(0, h.pos[indices[0]])
        self.zm.setPosition(1, h.pos[indices[1]])
        self.zm.setPosition(2, h.pos[indices[2]])
        for i in range(3, self.num_atoms):
            j = self.zm.entries[i].bond_index
            k = self.zm.entries[i].angle_index
            if k > 0:
                n = self.zm.entries[k].bond_index
            else:
                # First atom doesn't have a bond_index -- choose another atom for
                # torsion angle
                n = 0
                while n == j or n == k:
                    n += 1
            if -3 == h.atom_types[indices[i]].num_bonds and -3 != h.atom_types[indices[n]].num_bonds:
                # Ring atom: choose dihedral also on ring, if possible
                for q in range(i - 1, -1, -1):
                    if -3 == h.atom_types[indices[q]].num_bonds and (j != q and k != q) and (
                            self.zm.entries[q].bond_index == j or self.zm.entries[q].bond_index == k):
                        n = q
                        break
            self.zm.entries[i].dihedral_index = n
            v1 = h.pos[indices[i]] - h.pos[indices[j]]
            v2 = h.pos[indices[k]] - h.pos[indices[j]]
            n1 = np.cross(v1, v2)
            v1 = h.pos[indices[n]] - h.pos[indices[j]]
            n2 = np.cross(v1, v2)
            ang = RAD2DEG(calculate_angle(-n1, np.zeros(3), -n2))
            if np.isnan(ang):
                raise ValueError("Torsion angle for atom %d is NaN" % (i + 1))

        #define ANGLE_TARGET 10.0
        #ang = round(ang/ANGLE_TARGET)*ANGLE_TARGET;

            self.zm.entries[i].torsion_angle = ang
            zei = self.zm.entries[i]

            # Torsion corrections
            d2min = 1e9
            for j in range(360):
                v1 = self.zm.getPosition(i, v1)
                d1 = np.linalg.norm(v1 - h.pos[indices[i]])**2
                if d1 < TARGET_SQ_DIST:
                    break
                #printf("atom %d: sq_dist = %g\n", i+1, d1);
                if d1 < d2min:
                    d2min = d1
                    ang = self.zm.entries[i].torsion_angle
                self.zm.entries[i].torsion_angle += 1.0
                if self.zm.entries[i].torsion_angle > 360.0:
                    self.zm.entries[i].torsion_angle -= 360.0
            if 359 == j:
                self.zm.entries[i].torsion_angle = ang
            #if 0
            # for (j = 0; j < i; j++)
            # {
            #    zej = &m->zm->entries[j];
            #    if ((zei->bond_index == zej->bond_index) &&
            #        (zei->angle_index == zej->angle_index) &&
            #        (zei->dihedral_index == zej->dihedral_index) &&
            #        (fabs(adjustedAngle(zei->torsion_angle)-
            #              adjustedAngle(zej->torsion_angle)) < ANGLE_DIFF))
            #    {
            #       printf("adjusting torsion on atom %d from %g to ", i+1,
            #              zei->torsion_angle);

            #       zei->torsion_angle *= -1.0;
            #       zei->torsion_angle += ANGLE_DIFF;

            #       zei->torsion_angle = zej->torsion_angle + ANGLE_DIFF;
            #       printf("%g due to atom %d (torsion %g)\n", zei->torsion_angle,
            #              j+1, zej->torsion_angle);
            #       //break;
            #    }
            # }

            # if (zei->torsion_angle < 0.0)
            #    zei->torsion_angle = 360.0 + zei->torsion_angle;

            #endif
            #if 0
            # Now, some torsion angles may still be wrong, because acos() (in
            # getAngle()) returns the princpal value of the angle between Vectors;
            # this MAY NOT be a problem, if the new torsion angle does not put the
            # atom too close to previous atoms in the Monomer; if it does put the
            # atom too close, try 360-angle
            # getPosition(m->zm, i, &v1);
            # d1 = getElementR0(m->zm->entries[i].type.element_index);
            # for (j = 0; j < i; j++)
            # {
            #    getPosition(m->zm, j, &v2);
            #    d2 = getElementR0(m->zm->entries[j].type.element_index);
            #    if (getDist(&v1, &v2) < (d1+d2))
            #    {
            #       zei->torsion_angle = 360.0 - zei->torsion_angle;
            #       break;
            #    }
            # }
            # Adjust torsions close to 0 and 180 to enforce planar geometries
            # if ((zei->torsion_angle < 10.5) || (zei->torsion_angle > 349.5))
            #    zei->torsion_angle = 0.0;
            # if ((zei->torsion_angle < 190.5) && (zei->torsion_angle > 169.5))
            #    zei->torsion_angle = 180.0;
            #endif

        # Update bonds to reflect (possibly) new indices
        b = blist
        while b and hasattr(b, 'next'):
            b.index1 = rev_indices[b.index1]
            b.index2 = rev_indices[b.index2]
            b = b.next

        del indices
        del rev_indices
        del not_bb

    # ============================================================================
    # findExtraBonds()
    # ----------------------------------------------------------------------------
    # Result: for any bond in blist that is not found in m->zm, create a new Bond
    # and add it to m->extra_bonds and increment m->num_extra_bonds
    # ============================================================================
    def findExtraBonds(self, blist):
        b = blist

        while b and hasattr(b, 'next'):
            i1 = self.zm.entries[b.index1].bond_index
            i2 = self.zm.entries[b.index2].bond_index
            if b.index1 != i2 and b.index2 != i1:
                bnew = createBond(b.index1, b.index2)
                bnew.next = self.extra_bonds
                self.extra_bonds = bnew
                self.num_extra_bonds += 1
            b = b.next

    # ============================================================================
    # setFixedTorsion()
    # ----------------------------------------------------------------------------
    # Result: set the fixed torsion value for the indexth backbone atom of m
    # ============================================================================
    def setFixedTorsion(self, index, angle):
        self.torsion_angles[index] = {}
        self.torsion_angles[index][0] = angle

    # ============================================================================
    # readTorsionEnergies()
    # ----------------------------------------------------------------------------
    # Result: allocate and populate m->torsion_angles[index] and
    # m->torsion_probs[index] using the contents of the indicated file;
    # calls choke() on error
    # ============================================================================
    def readTorsionEnergies(self, index, path, temperature):
        s = createScanner(path)
        kT = kB * temperature

        # Expected format: "angle (deg)    energy (kcal/mol) \n"
        n = 0
        while TOK_EOF != s.getToken():
            n += 1
        n /= 2
        self.torsion_angles[index] = {}
        self.torsion_probs[index] = {}
        self.torsion_energies[index] = {}
        self.num_torsions[index] = n
        s.resetScanner()
        Emin = REAL_MAX
        for i in range(n):
            self.torsion_angles[index][i] = s.getRealToken()
            self.torsion_energies[index][i] = s.getRealToken()
            if self.torsion_energies[index][i] < Emin:
                Emin = self.torsion_energies[index][i]
        del s
        # Convert energies to Boltzmann probabilities, using temperature
        psum = 0.0
        for i in range(n):
            self.torsion_probs[index][i] = np.exp((Emin - self.torsion_energies[index][i]) / kT)
            psum += self.torsion_probs[index][i]
        # Normalize
        if psum > 1.0:
            for i in range(n):
                self.torsion_probs[index][i] /= psum
        # Store minimum probability (maximum energy) for the indexth torsion
        self.torsion_prob_min[index] = REAL_MAX
        for i in range(n):
            if self.torsion_probs[index][i] < self.torsion_prob_min[index]:
                self.torsion_prob_min[index] = self.torsion_probs[index][i]


# ============================================================================
# readMonomer()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Monomer populated
# by reading the indicated file (which may be .pdb, .zm, or .xyz); add
# AtomTypes; calls choke() on path errors (unknown extension, can't open file,
# etc).
# ============================================================================
def readMonomer(name, path, head_index, tail_index, eq_bond_scale):
    zm = {}
    blist = Bond()
    s = Scanner()

    # printf("\nMonomer %s...\n", name); fflush(stdout);
    # Read atomic species and positions into h
    p = path.rfind('.')
    if not p:
        raise NameError("Unable to determine file type of monomer file %s" % path)
    p += 1
    if path[p:] == "xyz":
        s = createScanner(path)
        h = readXYZ(s)
    elif path[p:] == "zm":
        if 0 != head_index:
            raise ValueError("For .zm input, head atom must be first")
        s = createScanner(path)
        h = readZM(s, zm)
    elif path[p:] == "pdb":
        h = readPDB(path)
    else:
        raise TypeError("Unable to read monomer file of type %s" % path[p:])
    if s:
        del s

    # Create Monomer
    m = Monomer()
    m.name = name
    m.num_atoms = h.num_atoms
    m.selection_count = 0
    m.create()

    # Fill m->bonds list and set m->num_bonds; each Bond's bond_type will be
    # set later...
    blist = findBonds(h, eq_bond_scale)

    # Determine atom types: element AND connectivity
    setAtomTypes(h, blist)

    # Now that we know the bonds and the AtomTypes, we can register the bond
    # types of each Bond in blist
    setBondTypes(h, blist)

    # Calculate internal coordinates, identify backbone atoms
    if zm:
        m.zm = zm
        m.num_bb = tail_index + 1
        for i in range(m.num_atoms):
            m.zm.entries[i].type = h.atom_types[i]
    else:
        m.zm = createZMatrix(m.num_atoms, 0)  # don't store positions
        m.fillZMatrix(h, blist, head_index, tail_index)

    # Identify any bonds in blist that are NOT represented in m->zm
    m.extra_bonds = Bond()
    m.num_extra_bonds = 0
    m.findExtraBonds(blist)
    del blist
    blist = Bond()

    # Setup torsions
    for i in range(m.num_bb):
        m.torsions[i] = Torsion(0)
        m.torsion_angles[i] = {}
        m.torsion_probs[i] = {}
        m.torsion_energies[i] = {}
        m.torsion_prob_min[i] = 0.0
        m.num_torsions[i] = 0

    # Masses
    m.head_mass = getElementMass(m.zm.entries[0].type.element_index)
    m.tail_mass = getElementMass(m.zm.entries[m.num_bb - 1].type.element_index)
    m.central_mass = 0.0
    for i in range(1, m.num_atoms):
        if (m.num_bb - 1) != i:
            m.central_mass += getElementMass(m.zm.entries[i].type.element_index)

    del h
    return m
