# ============================================================================
# system.py -- System functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

from .chain import Chain
from .domain import Domain, OccAtom
from .params import Params


class PolymerSystem:
    def __init__(self):
        self.chains = {}
        self.domains = {}
        self.rngs = {}
        self.pending_atoms = OccAtom()  # list of atoms that crossed Domain boundary
        self.pending_atoms.create()

    # ============================================================================
    # cleanupSystem()
    # ----------------------------------------------------------------------------
    # Result: free all fields of s
    # ============================================================================
    def cleanupSystem(self, p):
        i = 0

        if self.chains:
            for i in range(p.num_chais):
                del self.chains[i]
                self.chains[i] = []
            del self.chains
            self.chains = {}

        if self.domains:
            for i in range(p.total_domains):
                del self.domains[i]
                self.domains = {}
            del self.domains
            self.domains = {}

        if self.rngs:
            for i in range(p.total_domains):
                del self.rngs[i]
            del self.rngs
            self.rngs = {}

        oa1 = self.pending_atoms
        while oa1:
            oa2 = oa1
            del oa1
            oa1 = oa2
        self.pending_atoms = OccAtom()
        self.pending_atoms.create()
