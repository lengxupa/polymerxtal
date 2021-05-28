# ============================================================================
# vector.py -- Vector functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

# import numpy as np


# ============================================================================
# getNearestSqDist()
# ----------------------------------------------------------------------------
# Result: return the square distance between the nearest periodic images of v1
# and v2, given the dimensions (and half dimensions) of the system
# ============================================================================
def getNearestSqDist(v1, v2, sys_size, sys_half):

    diff = v1 - v2
    if diff[0] > sys_half[0]:
        diff[0] -= sys_size[0]
    if diff[0] < -sys_half[0]:
        diff[0] += sys_size[0]
    if diff[1] > sys_half[1]:
        diff[1] -= sys_size[1]
    if diff[1] < -sys_half[1]:
        diff[1] += sys_size[1]
    if diff[2] > sys_half[2]:
        diff[2] -= sys_size[2]
    if diff[2] < -sys_half[2]:
        diff[2] += sys_size[2]
    return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]
