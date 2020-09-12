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
from .utils import FREE


class Stereo:
    def __init__(self):
        self.name = ''
        self.pattern = 0  # flag: pattern or weight
        self.pattern_index = -1  # running index for a pattern
        self.curr_monomer = -1  # running index for added Monomers
        self.num_monomers = 0  # size of monomers[], weights[]
        self.selection_count = 0
        self.monomers = {}
        self.term = Monomer()
        self.term.create()
        self.weights = []

    def create(self):
        self.next = Stereo()
