"""
Functions to generate polymer crystals.
"""
import numpy as np

from polymerxtal.io import open_pdb

from .chain import Chain
from .molecule import calculate_density
from .move import Sphere, Cluster, Translator, Rotator


class Cell:
    def __init__(self, chain=Chain()):
        if not chain.built:
            chain.build_chain()
        self.chain = chain
        self.x = chain.x
        self.y = chain.y
        self.z = chain.z
        self.alpha = chain.alpha
        self.beta = chain.beta
        self.gamma = chain.gamma

    def hexagonal_packing(packing_density):
        self.crystal_name = self.chain.helix_name + "_hexagonal"

        symbols, coords = open_pdb(self.chain.helix_name + ".pdb")
        volume = self.x * self.y * self.z
        density = calculate_density(symbols, volume)

        L = np.sqrt(density * self.x * self.y / packing_density / np.sqrt(3))
        self.x = L * np.sqrt(3) / 2
        self.y = L / 2
