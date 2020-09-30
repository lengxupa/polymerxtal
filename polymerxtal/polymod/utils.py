# ============================================================================
# utils.py -- Utility functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import errno, os, random
import numpy as np

from .stdio import FILE


# ============================================================================
# openFile()
# ----------------------------------------------------------------------------
# Result: return a pointer to newly opened FILE; calls choke() to exit on
# failure
# ============================================================================
def openFile(path, mode):
    f = FILE()
    f.path = path
    f.file = open(path, mode)

    if not f.file:
        raise IOError("Unable to open %s: %s" % (path, os.strerror(errno.EIO)))
    return f


# ============================================================================
# selectWeight()
# ----------------------------------------------------------------------------
# Result: return a weight index in [0, num_weights) chosen according to
# weights, which sum to 1.0
# ============================================================================
def selectWeight(weights, num_weights, rng):
    i = 0
    target = random.random()
    p = 0.0

    while i < num_weights:
        if target >= p and target < (p + weights[i]):
            break
        p += weights[i]
        i += 1
    return i


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


# ============================================================================
# hashBin()
# ----------------------------------------------------------------------------
# Result: return the bin that holds pos, given the minimum position min, the
# size of the bins, and the number of bins in x and y
# ============================================================================
def hashBin(pos, mini, bin_size, nx, ny):
    cx = int((pos[0] - mini[0]) / bin_size[0])
    cy = int((pos[1] - mini[1]) / bin_size[1])
    cz = int((pos[2] - mini[2]) / bin_size[2])

    # Increase bin index first along x, then along y, then z
    return cx + nx * cy + nx * ny * cz


MIN_X_EDGE = 0
MIN_Y_EDGE = 1
MIN_Z_EDGE = 2
MAX_X_EDGE = 3
MAX_Y_EDGE = 4
MAX_Z_EDGE = 5


# ============================================================================
# getNeighborIndices()
# ----------------------------------------------------------------------------
# Result: fill indices, an int array of size 27, with the indices of regions
# that are neighbors to index in a system with Nx, Ny, Nz regions along each
# axis; set indices[i] to -1 if it is a repeated index; set edge_flags, an
# array of size 6, if index is on any edge
# ============================================================================
def getNeighborIndices(index, indices, edge_flags, Nx, Ny, Nz):
    Nxy = Nx * Ny
    gz = index / Nxy
    gy = (index - gz * Nxy) / Nx
    gx = index % Nx
    z_shift = Nz * Nxy
    i_shift = [[0, 0, 0], [-1, -1, 0], [0, -1, 0], [1, -1, 0], [-1, 0, 0], [1, 0, 0], [-1, 1, 0], [0, 1, 0], [1, 1, 0],
               [0, 0, -1], [-1, -1, -1], [0, -1, -1], [1, -1, -1], [-1, 0, -1], [1, 0, -1], [-1, 1, -1], [0, 1, -1],
               [1, 1, -1], [0, 0, 1], [-1, -1, 1], [0, -1, 1], [1, -1, 1], [-1, 0, 1], [1, 0, 1], [-1, 1, 1],
               [0, 1, 1], [1, 1, 1]]

    for i in range(6):
        edge_flags[i] = 0
    for i in range(27):
        indices[i] = index + i_shift[i][0] + i_shift[i][1] * Nx + i_shift[i][2] * Nxy
    if 0 == gx:
        # -X edge
        edge_flags[MIN_X_EDGE] += 1
        indices[1] += Nx
        indices[4] += Nx
        indices[6] += Nx
        indices[10] += Nx
        indices[13] += Nx
        indices[15] += Nx
        indices[19] += Nx
        indices[22] += Nx
        indices[24] += Nx
    if (Nx - 1) == gx:
        # +X edge
        edge_flags[MAX_X_EDGE] += 1
        indices[3] -= Nx
        indices[5] -= Nx
        indices[8] -= Nx
        indices[12] -= Nx
        indices[14] -= Nx
        indices[17] -= Nx
        indices[21] -= Nx
        indices[23] -= Nx
        indices[26] -= Nx
    if 0 == gy:
        # -Y edge
        edge_flags[MIN_Y_EDGE] += 1
        indices[1] += Nxy
        indices[2] += Nxy
        indices[3] += Nxy
        indices[10] += Nxy
        indices[11] += Nxy
        indices[12] += Nxy
        indices[19] += Nxy
        indices[20] += Nxy
        indices[21] += Nxy
    if (Ny - 1) == gy:
        # +Y edge
        edge_flags[MAX_Y_EDGE] += 1
        indices[6] -= Nxy
        indices[7] -= Nxy
        indices[8] -= Nxy
        indices[15] -= Nxy
        indices[16] -= Nxy
        indices[17] -= Nxy
        indices[24] -= Nxy
        indices[25] -= Nxy
        indices[26] -= Nxy
    if 0 == gz:
        # -Z edge
        edge_flags[MIN_Z_EDGE] += 1
        for i in range(9, 18):
            indices[i] += z_shift
    if (Nz - 1) == gz:
        # +Z edge
        edge_flags[MAX_Z_EDGE] += 1
        for i in range(18, 27):
            indices[i] -= z_shift

    # for (i = 0; i < 27; i++)
    # {
    #    for (j = i+1; j < 27; j++)
    #    {
    #       if ((indices[i] > -1) && (indices[i] == indices[j]))
    #          indices[j] = -1;
    #    }
    # }


# Conversions between radians and degrees
def RAD2DEG(rad):
    return rad * 57.295779513082325


def DEG2RAD(deg):
    return deg * 0.017453292519943


# Conversion from amu to grams
def AMU2GRAM(a):
    return a * 1.6605e-24


# Conversions between cubic Angstroms and cubic centimeters
def CA2CC(ca):
    return ca * 1.0e-24


def CC2CA(cc):
    return cc * 1.0e24


# Boltzmann constant in kcal/mol/K
kB = 0.0019872041
