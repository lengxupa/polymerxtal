# ============================================================================
# types.py -- AtomType typedef, prototypes for functions used in latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .element import getElementName
from .stdio import fputc  # FILE,


class AtomType:
    def __init__(self):
        self.element_index = -1
        self.num_bonds = -1  # negative value: resonant (-3), H-bond (-1)

    def __str__(self):
        return "element_index: %d\nnum_bonds: %d" % (self.element_index, self.num_bonds)

    # ============================================================================
    # writeAtomType()
    # ----------------------------------------------------------------------------
    # Result: write an AtomType string to the output FILE f
    # ============================================================================
    def writeAtomType(self, f):
        n = abs(self.num_bonds)

        f.printf(
            "%s with %d bond%s"
            % (getElementName(self.element_index), n, "s" if n > 1 else "")
        )
        if -3 == self.num_bonds:
            f.printf(" (resonant)")
        elif -1 == self.num_bonds:
            f.printf(" (H-bond)")

    # ============================================================================
    # writeAtomTypeDreiding()
    # ----------------------------------------------------------------------------
    # Result: write an AtomType Drediing string to the output FILE f
    # ============================================================================
    def writeAtomTypeDreiding(self, f):
        s = getElementName(self.element_index)
        NUM_FEL = 9
        fel = ["H", "F", "Cl", "Br", "I", "Na", "Ca", "Fe", "Zn"]

        f.printf("%s" % s)
        fputc("_", f)
        # ifdef HYDROGEN_BOND
        if -1 == self.num_bonds:
            f.printf("_HB")
        # endif
        i = 0
        while i < NUM_FEL:
            if s == fel[i]:
                break
            i += 1
        # If s is an entry in fel[], then we're done -- no hybridization needed
        if i == NUM_FEL:
            if self.num_bonds < 0:
                fputc("R", f)
            elif s == "O":
                f.printf("%d" % (self.num_bonds + 1))
            else:
                f.printf("%d" % (self.num_bonds - 1))


# Used only in this file
class AtomTypeList:
    def __init__(self):
        self.at = AtomType()
        self.index = -1

    def create(self):
        self.next = AtomTypeList()


def AT_EQUAL(at1, at2):
    return at1.element_index == at2.element_index and at1.num_bonds == at2.num_bonds


# Used only in this file
class BondTypeList:
    def __init__(self):
        self.at1 = AtomType()
        self.at2 = AtomType()
        self.index = -1

    def create(self):
        self.next = BondTypeList()


# File scope */
atl_top = AtomTypeList()
btl_top = BondTypeList()


# ============================================================================
# getAtomTypeIndex()
# ----------------------------------------------------------------------------
# Result: return the matching AtomType index or -1 if no match
# ============================================================================
def getAtomTypeIndex(t):
    atl = atl_top

    while atl and hasattr(atl, "next"):
        if AT_EQUAL(atl.at, t):
            return atl.index
        else:
            atl = atl.next
    return -1


# ============================================================================
# addAtomType()
# ----------------------------------------------------------------------------
# Result: register a new AtomType
# ============================================================================
def addAtomType(t):
    atl = AtomTypeList()
    atl.create()

    atl.at = t
    global atl_top
    atl.next = atl_top
    if atl_top and hasattr(atl_top, "next"):
        atl.index = atl_top.index + 1
    else:
        atl.index = 1
    atl_top = atl


# ============================================================================
# writeAtomTypes()
# ----------------------------------------------------------------------------
# Result: write a formatted list of AtomTypes, element and number of bonds, to
# the output FILE f
# ============================================================================
def writeAtomTypes(f):
    atl = atl_top

    f.printf(
        "%d known atom types:\n"
        % (atl_top.index if atl_top and hasattr(atl_top, "next") else 0)
    )
    while atl and hasattr(atl, "next"):
        f.printf("   ")
        atl.at.writeAtomType(f)
        f.printf("\n")
        atl = atl.next


# ============================================================================
# writeAtomTypesDreiding()
# ----------------------------------------------------------------------------
# Result: write a list of AtomTypes, index and Dreiding symbol, to the output
# FILE f
# ============================================================================
def writeAtomTypesDreiding(f):
    atl = atl_top

    while atl and hasattr(atl, "next"):
        f.printf("%d  " % atl.index)
        atl.at.writeAtomTypeDreiding(f)
        f.printf("\n")
        atl = atl.next


# ============================================================================
# getBondTypeIndex()
# ----------------------------------------------------------------------------
# Result: return the matching bond type index or -1 if no match
# ============================================================================
def getBondTypeIndex(at1, at2):
    btl = btl_top

    while btl and hasattr(btl, "next"):
        if (AT_EQUAL(at1, btl.at1) and AT_EQUAL(at2, btl.at2)) or (
            AT_EQUAL(at2, btl.at1) and AT_EQUAL(at1, btl.at2)
        ):
            return btl.index
        else:
            btl = btl.next
    return -1


# ============================================================================
# addBondType()
# ----------------------------------------------------------------------------
# Result: register a new bond type
# ============================================================================
def addBondType(at1, at2):
    btl = BondTypeList()
    btl.create()

    btl.at1 = at1
    btl.at2 = at2
    global btl_top
    btl.next = btl_top
    if btl_top and hasattr(btl_top, "next"):
        btl.index = btl_top.index + 1
    else:
        btl.index = 1
    btl_top = btl


# ============================================================================
# writeBondTypes()
# ----------------------------------------------------------------------------
# Result: write a formatted list of bond types, between two AtomTypes, to the
# output FILE f
# ============================================================================
def writeBondTypes(f):
    btl = btl_top

    f.printf(
        "%d known bond types:\n"
        % (btl_top.index if (btl_top and hasattr(btl_top, "next")) else 0)
    )
    while btl and hasattr(btl, "next"):
        f.printf("   (")
        btl.at1.writeAtomType(f)
        f.printf(") -- (")
        btl.at2.writeAtomType(f)
        f.printf(")\n")
        btl = btl.next


# ============================================================================
# writeBondTypesDreiding()
# ----------------------------------------------------------------------------
# Result: write a list of BondTypes, index and Dreiding atom types, to the
# output FILE f
# ============================================================================
def writeBondTypesDreiding(f):
    btl = btl_top

    while btl and hasattr(btl, "next"):
        f.printf("%d  " % btl.index)
        btl.at1.writeAtomTypeDreiding(f)
        f.printf("  ")
        btl.at2.writeAtomTypeDreiding(f)
        f.printf("\n")
        btl = btl.next
