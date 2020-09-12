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


# Excluded cylindrical region
class ExclCylinder:
    def __init__(self):
        self.start = np.zeros(3)  # starting point in center of one end
        self.axis = np.zeros(3)
        self.radius = 0.  # Angstroms
        self.length = 0.  # Angstroms
        self.invert = 0

    # ============================================================================
    # createExclCylinder()
    # ----------------------------------------------------------------------------
    # Result: return a pointer to a newly allocated, initialized ExclCylinder; all
    # fields are filled later; call choke() on error
    # ============================================================================
    def create(self):
        self.invert = 0
        self.next = ExclCylinder()


# Excluded rectangular region
class ExclSlab:
    def __init__(self):
        self.min = np.zeros(3)
        self.max = np.zeros(3)
        self.invert = 0

    # ============================================================================
    # createExclSlab()
    # ----------------------------------------------------------------------------
    # Result: return a pointer to a newly allocated, initialized ExclSlab; all
    # fields are filled later; call choke() on error
    # ============================================================================
    def create(self):
        self.invert = 0
        self.next = ExclSlab()


# Excluded spherical region
class ExclSphere:
    def __init__(self):
        self.center = np.zeros(3)
        self.radius = 0.  # Angstroms
        self.invert = 0

    # ============================================================================
    # createExclSphere()
    # ----------------------------------------------------------------------------
    # Result: return a pointer to a newly allocated, initialized ExclSphere; all
    # fields are filled later; call choke() on error
    # ============================================================================
    def create(self):
        self.invert = 0
        self.next = ExclSphere()
