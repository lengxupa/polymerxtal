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

from .zmatrix import ZMatrix, Vector, createZMatrix
from .utils import FREE
from .types import AtomType
from .scan import *
from .stdio import rewind


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
        h.pos[i] = Vector()
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
    return eval(string)


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
        h.pos[i].x = s.getRealToken()
        h.pos[i].y = s.getRealToken()
        h.pos[i].z = s.getRealToken()
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
    v0 = Vector()

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
        zmptr.getPosition(i, h.pos[i])

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
    h = {}
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

    f = open(path, "r")
    n = 0
    line = f.readline()
    while line:
        if (line[:6] == "ATOM  ") or (line[:6] == "HETATM"):
            n += 1
    h = createHolder(n)

    rewind(f)
    i = 0
    line = f.readline()
    while line:
        leng = len(line)
        if (line[:6] == "ATOM  ") or (line[:6] == "HETATM"):
            h.pos[i].x = eval(line[30:38])
            h.pos[i].y = eval(line[38:46])
            h.pos[i].z = eval(line[46:54])
            # Check element symbol first...
            if (not line[77].isspace()) and leng > 77:
                line[78] = '\0'
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
    f.close()
    return h


# Entry in a linked list
class Monomer:
    def __init__(self):
        self.name = ''
        self.num_atoms = 0
        self.num_extra_bonds = 0
        self.num_bb = 0  # Number of backbone atoms; head: 0;  tail: num_bb-1
        self.selection_count = 0
        self.zm = []
        self.extra_bonds = Bond()
        self.extra_bonds.create()
        self.torsions = []
        self.torsion_angles = {}
        self.torsion_energies = {}
        self.torsion_probs = {}
        self.torsion_prob_min = []
        self.head_mass = 0.
        self.tail_mass = 0.
        self.central_mass = 0.  # all atoms except head and tail
        self.num_torsions = []

    def create(self):
        self.next = Monomer()


# ============================================================================
# readMonomer()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Monomer populated
# by reading the indicated file (which may be .pdb, .zm, or .xyz); add
# AtomTypes; calls choke() on path errors (unknown extension, can't open file,
# etc).
# ============================================================================
def readMonomer(name, path, head_index, tail_index, eq_bond_scale):
    m = Monomer()
    m.create()
    h = {}
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
        s = s.createScanner(path)
        h = readXYZ(s)
    elif path[p:] == "zm":
        if 0 != head_index:
            raise ValueError("For .zm input, head atom must be first")
        s = s.createScanner(path)
        h = readZM(s, zm)
    elif path[p:] == "pdb":
        h = readPDB(path)
    else:
        raise TypeError("Unable to read monomer file of type %s" % path[p:])
    if s:
        del s

    # Unfinished


#Monomer *
#readMonomer(const char * restrict name, const char * restrict path,
#            int head_index, int tail_index, Real eq_bond_scale)
#{
#   Monomer *m;
#   Holder *h;
#   ZMatrix *zm = NULL;
#   Bond *blist = NULL;
#   Scanner *s = NULL;
#   char *p;
#   int i;


# ============================================================================
# createBond()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Bond; calls
# choke() if allocation fails
# ============================================================================
def createBond(index1, index2):
    b = Bond()
    b.create()

    b.index1 = index1
    b.index2 = index2
    b.visited = 0
    b.next = Bond()
    return b
