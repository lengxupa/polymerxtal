# ============================================================================
# utils.py -- Utility functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .vector import Vector


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
    while pos.x < system_min.x:
        pos.x += system_size.x
    while pos.x > system_max.x:
        pos.x -= system_size.x
    if pos.x == system_max.x:
        pos.x = system_min.x

    while pos.y < system_min.y:
        pos.y += system_size.y
    while pos.y > system_max.y:
        pos.y -= system_size.y
    if pos.y == system_max.y:
        pos.y = system_min.y

    while pos.z < system_min.z:
        pos.z += system_size.z
    while pos.z > system_max.z:
        pos.z += system_size.z
    if pos.z == system_max.z:
        pos.z = system_min.z


# Conversions between radians and degrees
def DEG2RAD(deg):
    return deg * 0.017453292519943
