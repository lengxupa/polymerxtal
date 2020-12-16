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

from .params import Params
from .utils import getNeighborIndices, hashBin


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
        self.nbr_bins = {}
        self.domain_edge = {}
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
        self.nbr_bins = {}  # Used for all Bins if not storing neighbors
        self.domain_edge = {}  # Used for all Bins if not storing neighbors
        self.nbr_domains = {}  # always stored: relatively few domains
        self.system_edge = {}  # always stored: relatively few domains
        self.bins = {}

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

    # ============================================================================
    # addAtom()
    # ----------------------------------------------------------------------------
    # Result: add a new OccAtom to the Domain d
    # ============================================================================
    def addAtom(self, pos, chain, atom):
        oa = OccAtom()
        oa.create()
        n = hashBin(pos, self.min, self.bin_size, self.num_bins_x, self.num_bins_y)

        oa.chain = chain
        oa.atom = atom
        oa.next = self.bins[n].oa_list
        self.bins[n].oa_list = oa

    # ============================================================================
    # addOccAtom()
    # ----------------------------------------------------------------------------
    # Result: add an existing OccAtom to the Domain d
    # ============================================================================
    def addOccAtom(self, pos, oa):
        n = hashBin(pos, self.min, self.bin_size, self.num_bins_x, self.num_bins_y)

        oa.next = self.bins[n].oa_list
        self.bins[n].oa_list = oa

    # ============================================================================
    # removeDeadChains()
    # ----------------------------------------------------------------------------
    # Result: remove all OccAtoms with matching chain_index and atom index >=
    # atom_start from the Domain d
    # ============================================================================
    def removeDeadChains(self, chain_index, atom_start):

        for i in range(self.total_bins):
            oa = self.bins[i].oa_list
            prev = OccAtom()
            while oa and hasattr(os, "next"):
                if oa.chain == chain_index and oa.atom >= atom_start:
                    if prev and hasattr(prev, "next"):
                        prev.next = oa.next
                    else:  # top of list
                        self.bins[i].oa_list = oa.next
                    oa.next = OccAtom()
                prev = oa
                oa = oa.next


# ============================================================================
# createDomain()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Domain
# ============================================================================
def createDomain(index, mini, maxi, p):
    d = Domain()

    # printf("Creating domain %d from (%g, %g, %g) to (%g, %g, %g)\n", index,
    # min->x, min->y, min->z, max->x, max->y, max->z);
    # fflush(stdout);

    d.index = index
    d.min = mini
    d.max = maxi
    # Bin size is the same in all dimensions; make it a Vector in order to
    # use hashBin() (hashBin() also hashes Domains, which may not be same
    # in each dimension)

    d.bin_size[0] = p.grid_size
    d.bin_size[1] = p.grid_size
    d.bin_size[2] = p.grid_size
    d.num_bins_x = int(np.ceil((maxi[0] - mini[0]) / p.grid_size))
    d.num_bins_y = int(np.ceil((maxi[1] - mini[1]) / p.grid_size))
    d.num_bins_z = int(np.ceil((maxi[2] - mini[2]) / p.grid_size))
    d.total_bins = d.num_bins_x * d.num_bins_y * d.num_bins_z
    for i in range(d.total_bins):
        d.bins[i] = Bin()
    if p.recalculate_neighbors:
        # Fill these arrays for each Bin at each iteration
        for i in range(27):
            d.nbr_bins[i] = 0
        for i in range(6):
            d.domain_edge[i] = 0
    else:
        d.nbr_bins = {}
        d.domain_edge = {}
    for i in range(d.total_bins):
        if p.recalculate_neighbors:
            d.bins[i].nbr_bins = {}
            d.bins[i].domain_edge = {}
        else:
            for j in range(27):
                d.bins[i].nbr_bins[j] = 0
            for j in range(6):
                d.bins[i].domain_edge[j] = 0
            getNeighborIndices(
                i,
                d.bins[i].nbr_bins,
                d.bins[i].domain_edge,
                d.num_bins_x,
                d.num_bins_y,
                d.num_bins_z,
            )
        d.bins[i].oa_list = OccAtom()
    getNeighborIndices(
        index,
        d.nbr_domains,
        d.system_edge,
        p.num_domains_x,
        p.num_domains_y,
        p.num_domains_z,
    )
    return d
