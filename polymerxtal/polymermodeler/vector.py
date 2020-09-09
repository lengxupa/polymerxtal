# ============================================================================
# vector.py -- Vector typedef, prototypes for functions used in latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np


class Vector:
    def __init__(self):
        self.x = 0.
        self.y = 0.
        self.z = 0.

    # ============================================================================
    # addVectors()
    # ----------------------------------------------------------------------------
    # Result: set result = lhs + rhs
    # ============================================================================

    def addVectors(self, lhs, rhs):
        self.x = lhs.x + rhs.x
        self.y = lhs.y + rhs.y
        self.z = lhs.z + rhs.z

    # ============================================================================
    # subtractVectors()
    # ----------------------------------------------------------------------------
    # Result: set result = lhs - rhs
    # ============================================================================

    def subtractVectors(self, lhs, rhs):
        self.x = lhs.x - rhs.x
        self.y = lhs.y - rhs.y
        self.z = lhs.z - rhs.z

    # ============================================================================
    # scaleVector()
    # ----------------------------------------------------------------------------
    # Result: scale each component of v by scale
    # ============================================================================

    def scaleVector(self, scale):
        self.x *= scale
        self.y *= scale
        self.z *= scale

    # ============================================================================
    # crossVectors()
    # ----------------------------------------------------------------------------
    # Result: set result = lhs x rhs (vector product)
    # ============================================================================

    def crossVectors(self, lhs, rhs):
        self.x = lhs.y * rhs.z - lhs.z * rhs.y
        self.y = lhs.z * rhs.x - lhs.x * rhs.z
        self.z = lhs.x * rhs.y - lhs.y * rhs.x

    # ============================================================================
    # getLength()
    # ----------------------------------------------------------------------------
    # Result: return the length of v (assuming a Vector from the origin)
    # ============================================================================

    def getLength(self):
        return np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    # ============================================================================
    # normalizeVector()
    # ----------------------------------------------------------------------------
    # Result: adjust components of v so that |v| = 1
    # ============================================================================

    def normalizeVector(self):
        leng = 1.0 / self.getLength()

        self.x *= leng
        self.y *= leng
        self.z *= leng
