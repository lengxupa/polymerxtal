# ============================================================================
# system.py -- System functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

from .chain import Chain
from .domain import Domain, OccAtom
from .params import Params
from .utils import foldPosition, hashBin
from .zmatrix import ZMatrix


class PolymerSystem:
    # ============================================================================
    # initSystem()
    # ----------------------------------------------------------------------------
    # Result: initialize the fields of s to NULL
    # ============================================================================
    def __init__(self):
        self.chains = {}
        self.domains = {}
        self.rngs = {}
        self.pending_atoms = OccAtom()  # list of atoms that crossed Domain boundary

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

    # ============================================================================
    # addPendingAtom()
    # ----------------------------------------------------------------------------
    # Result: add a pending atom (an atom whose position is controlled by one
    # Domain but has crossed the boundary to another Domain)
    # ============================================================================
    def addPendingAtom(self, chain, atom):
        oa = OccAtom()
        oa.create()

        oa.chain = chain
        oa.atom = atom
        oa.next = self.pending_atoms
        self.pending_atoms = oa

    # ============================================================================
    # getPendingAtoms()
    # ----------------------------------------------------------------------------
    # Result: add pending atoms for the Domain d and remove from s->pending_atoms
    # ============================================================================
    def getPendingAtoms(self, d, p):
        oa = self.pending_atoms
        prev = OccAtom()

        while oa and hasattr(oa, 'next'):
            next = oa.next
            pos = self.chains[oa.chain].zm.getPosition(oa.atom)
            foldPosition(pos, p.system_min, p.system_max, p.system_size)
            domain = hashBin(pos, p.system_min, p.domain_size, p.num_domains_x, p.num_domains_y)
            if domain == d.index:
                if prev and hasattr(prev, 'next'):
                    prev.next = next
                else:
                    self.pending_atoms = next
                d.addOccAtom(pos, oa)
            else:
                prev = oa
            oa = next
