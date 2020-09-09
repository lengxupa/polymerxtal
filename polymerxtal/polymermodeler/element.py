# ============================================================================
# element.py -- Element functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .scan import Scanner, TOK_EOF

# File scope

elements = {}
num_elements = 0


class Element:
    def __init__(self):
        self.name = ''
        self.mass = 0.  # amu
        self.R0 = 0.  # equilibrium bond length; A
        self.LJ_R0 = 0.  # Lennard-Jones length parameter; A
        self.LJ_D0 = 0.  # Lennard-Jones energy parameter; kcal/mol


# ============================================================================
# readElements()
# ----------------------------------------------------------------------------
# Result: read element data from the named file; calls choke() on error
# ============================================================================


def readElements(path):
    s = Scanner()
    s.createScanner(path)
    num_tokens = 0

    while TOK_EOF != s.getToken():
        num_tokens += 1
    s.resetScanner()
    num_elements = num_tokens / 5
    for n in range(num_elements):
        elements[n] = Element()
        s.getToken()
        elements[n].name = s.tokstr
        elements[n].mass = s.getRealToken()
        elements[n].R0 = s.getRealToken()
        elements[n].LJ_R0 = s.getRealToken()
        elements[n].LJ_D0 = s.getRealToken()
    del s


# ============================================================================
# getElement()
# ----------------------------------------------------------------------------
# Result: return an Element pointer for the specified index; calls choke() if
# index is out of range
# ============================================================================


def getElement(index):
    if (index < 0) or (index >= num_elements):
        raise ValueError("Element index %d is out of range" % index)
    return elements[index]


# ============================================================================
# getElementName()
# ----------------------------------------------------------------------------
# Result: return a const char pointer to the specified Element name; calls
# choke() if index is out of range
# ============================================================================


def getElementName(index):
    el = getElement(index)

    return el.name


# ============================================================================
# getElementLJ_R0()
# ----------------------------------------------------------------------------
# Result: return the Lennard-Jones length parameter for the specified Element;
# calls choke() if index is out of range
# ============================================================================


def getElementLJ_R0(index):
    el = getElement(index)

    return el.LJ_R0


# ============================================================================
# getElementLJ_D0()
# ----------------------------------------------------------------------------
# Result: return the Lennard-Jones energy parameter for the specified Element;
# calls choke() if index is out of range
# ============================================================================


def getElementLJ_D0(index):
    el = getElement(index)

    return el.LJ_D0
