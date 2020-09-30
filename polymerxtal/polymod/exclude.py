# ============================================================================
# exclude.py -- ExclCylinder and ExclSlab typedefs, prototypes for functions
# used in latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

from polymerxtal import calculate_angle


# Excluded cylindrical region
class ExclCylinder:
    def __init__(self):
        self.start = np.zeros(3)  # starting point in center of one end
        self.axis = np.zeros(3)
        self.radius = 0.  # Angstroms
        self.length = 0.  # Angstroms
        self.invert = 0

    def create(self):
        self.next = ExclCylinder()

    # ============================================================================
    # insideExclCylinder()
    # ----------------------------------------------------------------------------
    # Result: return 1 if the position Vector p is inside ec, else return 0
    # ============================================================================
    def insideExclCylinder(self, p):

        v = p - self.start
        return np.linalg.norm(p - self.start) * np.sin(calculate_angle(-v, np.zeros(3), -self.axis)) < self.radius


# ============================================================================
# createExclCylinder()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized ExclCylinder; all
# fields are filled later; call choke() on error
# ============================================================================
def createExclCylinder():
    ec = ExclCylinder()
    ec.create()

    ec.invert = 0
    ec.next = ExclCylinder()
    return ec


# Excluded rectangular region
class ExclSlab:
    def __init__(self):
        self.min = np.zeros(3)
        self.max = np.zeros(3)
        self.invert = 0

    def create(self):
        self.next = ExclSlab()

    # ============================================================================
    # insideExclSlab()
    # ----------------------------------------------------------------------------
    # Result: return 1 if the position Vector p is inside es, else return 0
    # ============================================================================
    def insideExclSlab(self, p):
        return p[0] > self.min[0] and p[0] < self.max[0] and p[1] > self.min[1] and p[1] < self.max[1] and p[
            2] > self.min[2] and p[2] < self.max[2]


# ============================================================================
# createExclSlab()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized ExclSlab; all
# fields are filled later; call choke() on error
# ============================================================================
def createExclSlab():
    es = ExclSlab()
    es.create()

    es.invert = 0
    es.next = ExclSlab()
    return es


# Excluded spherical region
class ExclSphere:
    def __init__(self):
        self.center = np.zeros(3)
        self.radius = 0.  # Angstroms
        self.invert = 0

    def create(self):
        self.next = ExclSphere()

    # ============================================================================
    # insideExclSphere()
    # ----------------------------------------------------------------------------
    # Result: return 1 if the position Vector p is inside es, else return 0
    # ============================================================================
    def insideExclSphere(self, p):
        return np.linalg.norm(p - self.center)**2 < self.radius * self.radius


# ============================================================================
# createExclSphere()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized ExclSphere; all
# fields are filled later; call choke() on error
# ============================================================================
def createExclSphere():
    es = ExclSphere()
    es.create()

    es.invert = 0
    es.next = ExclSphere()
    return es
