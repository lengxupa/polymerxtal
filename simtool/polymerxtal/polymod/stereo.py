# ============================================================================
# stereo.py -- Stereo typedef, prototypes for functions used in latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .monomer import Monomer
from .utils import selectWeight


class Stereo:
    def __init__(self):
        self.name = ""
        self.pattern = 0  # flag: pattern or weight
        self.pattern_index = -1  # running index for a pattern
        self.curr_monomer = -1  # running index for added Monomers
        self.num_monomers = 0  # size of monomers[], weights[]
        self.selection_count = 0
        self.monomers = {}
        self.term = Monomer()
        self.term.create()
        self.weights = {}

    def create(self):
        self.next = Stereo()

    # ============================================================================
    # addStereoMonomer()
    # ----------------------------------------------------------------------------
    # Result: add a new Monomer with indicated weight; calls choke() if the
    # number of Monomers passed to createStereo() has been reached already
    # ============================================================================
    def addStereoMonomer(self, m, weight):
        if self.curr_monomer < self.num_monomers:
            self.monomers[self.curr_monomer] = m
            self.weights[self.curr_monomer] = weight
            self.curr_monomer += 1
        else:
            raise ValueError("Invalid Stereo monomer index %d" % self.curr_monomer)

    # ============================================================================
    # getNextMonomer()
    # ----------------------------------------------------------------------------
    # Result: return a pointer to the next Monomer to be added to a Chain
    # ============================================================================
    def getNextMonomer(self, rng):

        if self.pattern:
            m = self.monomers[self.pattern_index]
            self.pattern_index += 1
            if self.pattern_index == self.curr_monomer:
                self.pattern_index = 0
        else:
            m = self.monomers[selectWeight(self.weights, self.num_monomers, rng)]
        return m


# ============================================================================
# createStereo()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Stereo
# (monomers[] and weights[] are filled later); calls choke() on allocation
# error
# ============================================================================
def createStereo(name, pattern, num_monomers):
    s = Stereo()
    s.create()

    s.name = name
    s.pattern = pattern
    s.pattern_index = 0  # increment with calls to getNextMonomer()
    s.curr_monomer = 0  # increment with calls to addStereoMonomer()
    s.num_monomers = num_monomers
    for i in range(num_monomers):
        s.monomers[i] = Monomer()
        s.weights[i] = 0.0
    s.term = Monomer()
    s.next = Stereo()
    s.selection_count = 0
    return s
