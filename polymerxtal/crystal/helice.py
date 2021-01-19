"""
This module is for functions Helice class.
"""


class Helice:
    def __init__(self, atoms=2, motifs=1, turns=1, sub_type=0):
        self.atoms = atoms
        self.motifs = motifs
        self.turns = turns
        self.sub_type = sub_type

    def __str__(self):
        return "%d_%d_%d" % (self.atoms, self.motifs, self.turns)
