# ============================================================================
# element.py -- Element functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .scan import Scanner, createScanner, TOK_EOF

# File scope
elements = {}
num_elements = 0


class Element:
    def __init__(self):
        self.name = ""
        self.mass = 0.0  # amu
        self.R0 = 0.0  # equilibrium bond length; A
        self.LJ_R0 = 0.0  # Lennard-Jones length parameter; A
        self.LJ_D0 = 0.0  # Lennard-Jones energy parameter; kcal/mol


# ============================================================================
# readElements()
# ----------------------------------------------------------------------------
# Result: read element data from the named file; calls choke() on error
# ============================================================================
def readElements(path):
    s = createScanner(path)
    num_tokens = 0

    while TOK_EOF != s.getToken():
        num_tokens += 1
    s.resetScanner()
    global num_elements
    num_elements = int(num_tokens / 5)
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
# getElementIndex()
# ----------------------------------------------------------------------------
# Result: return the index of the named Element; calls choke() if no match
# ============================================================================
def getElementIndex(name):
    n = 0

    while n < num_elements:
        if name == elements[n].name:
            break
        n += 1
    if n == num_elements:
        raise NameError('Unknown element: "%s"' % name)
    return n


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
# getElementMass()
# ----------------------------------------------------------------------------
# Result: return the atomic mass for the specified Element; calls choke() if
# index is out of range
# ============================================================================
def getElementMass(index):
    el = getElement(index)

    return el.mass


# ============================================================================
# getElementR0()
# ----------------------------------------------------------------------------
# Result: return the equilibrium bond length for the specified Element; calls
# choke() if index is out of range
# ============================================================================
def getElementR0(index):
    el = getElement(index)

    return el.R0


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
