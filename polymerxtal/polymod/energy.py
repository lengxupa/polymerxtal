# ============================================================================
# energy.py -- Energy expressions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

from .element import *  # getElement*()

# File scope
self_avoid_cutoff_sq = 0.0


# ============================================================================
# setSelfAvoidCutoff()
# ----------------------------------------------------------------------------
# Result: set square self-avoiding cutoff energy
# ============================================================================
def setSelfAvoidCutoff(rcut):
    self_avoid_cutoff_sq = rcut * rcut


class Energy:

    # ============================================================================
    # energyLJ()
    # ----------------------------------------------------------------------------
    # Result: return the Lennard-Jones interaction energy for a pair of atoms of
    # types type1 and type2, separated by the square distance r2
    # ============================================================================
    def energyLJ(type1, type2, r2):
        p2 = getElementLJ_R0(type1) * getElementLJ_R0(type2) / r2  # (R/R0)^-2
        p6 = p2 * p2 * p2

        return np.sqrt(getElementLJ_D0(type1) * getElementLJ_D0(type2)) * p6 * (p6 - 2.0)

    # ============================================================================
    # energyHardCore()
    # ----------------------------------------------------------------------------
    # Result: return a large value if r2 is less the square equilibrium bond
    # distance for types type1 and type2, else return 0
    # ============================================================================
    def energyHardCore(type1, type2, r2):
        r = getElementR0(type1) + getElementR0(type2)

        return 1e3 if (r2 < r * r) else 0.0

    # ============================================================================
    # energySelfAvoid()
    # ----------------------------------------------------------------------------
    # Result: return a large value if r2 is less than the square cutoff, else
    # return 0
    # ============================================================================
    def energySelfAvoid(type1, type2, r2):
        return 1e3 if (r2 < self_avoid_cutoff_sq) else 0.0

    # ============================================================================
    # energyNone()
    # ----------------------------------------------------------------------------
    # Result: return 0.0 for all cases
    # ============================================================================
    def energyNone(type1, type2, r2):
        return 0.0
