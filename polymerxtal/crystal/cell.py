"""
Functions to generate polymer crystals.
"""
import numpy as np

from polymerxtal.io import open_pdb

from .chain import Chain
from .molecule import calculate_density
from .move import Sphere, Cluster, Translator, Rotator


class Cell:
    def __init__(
        self, **kwargs
    ):  # ,chain=Chain(), packing='Hexagonal', density=0.857):
        self.chain = Chain()
        if not chain.built:
            chain.build_chain()
        self.chain = chain
        self.a = chain.x
        self.b = chain.y
        self.c = chain.z
        self.alpha = chain.alpha
        self.beta = chain.beta
        self.gamma = chain.gamma

        self.packing = packing
        self.density = density

    def build_crystal(packing="Hexagonal", density=0.857):
        self.crystal_name = self.chain.helix_name + "_hexagonal"
        self.packing = packing
        self.density = density

        symbols, coords = open_pdb(self.chain.helix_name + ".pdb")
        volume = self.x * self.y * self.z
        den = calculate_density(symbols, volume)

        if packing == "Hexagonal":
            L = np.sqrt(den * self.x * self.y / density / np.sqrt(3))
            self.x = L * np.sqrt(3) / 2
            self.y = L / 2
        else:
            raise Exception(
                "Unknown packing method %s\nCurrently available method includes Hexagonal"
                % packing
            )
