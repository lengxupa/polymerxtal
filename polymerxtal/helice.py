"""
This module is for functions Helice class.
"""

class Helice:
    def __init__(self, atoms=2, motifs=3, turns=1):
        self.atoms = atoms
        self.motifs = motifs
        self.turns = turns

    def __str__(self):
        return f'{self.atoms}_{self.motifs}_{self.turns}'
