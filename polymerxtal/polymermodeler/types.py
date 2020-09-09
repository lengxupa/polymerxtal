# ============================================================================
# types.py -- AtomType typedef, prototypes for functions used in latch module
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================


class AtomType:
    def __init__(self):
        self.element_index = -1
        self.num_bonds = -1  # negative value: resonant (-3), H-bond (-1)

    def __str__(self):
        return f'element_index: {self.element_index}\nnum_bonds: {self.num_bonds}'
