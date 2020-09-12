# ============================================================================
# domain.py -- Domain functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

from .utils import FREE


# Spatial location of an atom: chain and atom indices
class OccAtom:
    def __init__(self):
        self.chain = 0
        self.atom = 0

    def create(self):
        self.next = OccAtom()


# Smallest spatial region
class Bin:
    def __init__(self):
        self.nbr_bins = []
        self.domain_edge = []
        self.oa_list = OccAtom()
        self.oa_list.create()


# Spatial region divided into bins: identify neighboring atoms
class Domain:
    def __init__(self):
        self.min = np.zeros(3)
        self.max = np.zeros(3)
        self.bin_size = np.zeros(3)
        self.index = -1
        self.num_bins_x = 0
        self.num_bins_y = 0
        self.num_bins_z = 0
        self.total_bins = 0
        self.nbr_bins = []  # Used for all Bins if not storing neighbors
        self.domain_edge = []  # Used for all Bins if not storing neighbors
        self.nbr_domains = []  # always stored: relatively few domains
        self.system_edge = []  # always stored: relatively few domains
        self.bins = []

    # ============================================================================
    # clearDomain()
    # ----------------------------------------------------------------------------
    # Result: remove all OccAtoms from the Domain d
    # ============================================================================
    def clearDomain(self):
        for i in range(self.total_bins):
            oa = self.bins[i].oa_list
            while oa:
                next = oa.next
                del oa.next
                del oa
                oa = next
            del self.bins[i].oa_list
