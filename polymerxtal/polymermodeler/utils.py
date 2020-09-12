# ============================================================================
# utils.py -- Utility functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np


# Free a list and set it to [] */
def FREE(ptr):
    if ptr:
        ptr = []


# ============================================================================
# foldPosition()
# ----------------------------------------------------------------------------
# Result: fold pos until it is inside (system_min, system_max)
# ============================================================================
def foldPosition(pos, system_min, system_max, system_size):
    while pos[0] < system_min[0]:
        pos[0] += system_size[0]
    while pos[0] > system_max[0]:
        pos[0] -= system_size[0]
    if pos[0] == system_max[0]:
        pos[0] = system_min[0]

    while pos[1] < system_min[1]:
        pos[1] += system_size[1]
    while pos[1] > system_max[1]:
        pos[1] -= system_size[1]
    if pos[1] == system_max[1]:
        pos[1] = system_min[1]

    while pos[2] < system_min[2]:
        pos[2] += system_size[2]
    while pos[2] > system_max[2]:
        pos[2] += system_size[2]
    if pos[2] == system_max[2]:
        pos[2] = system_min[2]


# Conversions between radians and degrees
def DEG2RAD(deg):
    return deg * 0.017453292519943
