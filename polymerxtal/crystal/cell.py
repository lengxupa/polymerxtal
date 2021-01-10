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
        """polymerxtal.crystal.cell.Cell class contains all the information about the unit cell of polymer crystal

        Parameters
        ----------
        polymer_type : str (optional)
            Current support polymer types include:
                PAN - Polyacrylonitrile
                PE - Polyethylene
                PP - Polypropylene
                PS - Polystyrene
                POM - Polyoxymethylene
                PTFE - Polytetrafluoroethylene
                PVC - Polyvinyl chloride
        helice : Helice class (optional)
            Current support helicities include:
                Helice(2,3,1)
                Helice(2,2,1)
                Helice(2,1,1) - Planar zig-zag
        num_monomers : int (optional)
            Number of monomers
        tacticity : str (optional)
            Tacticity of the polymer chain. Current support tacticities include: isotactic, atactic, syndiotactic
        chiriality : str (optional)
            Chiriality of the polymer chain. Current support chirialities include: left, right
        head_tail_defect_ratio : float (optional)
            Ratio of head-to-head and tail-to-tail connections
        configs : int (optional)
            Number of attempts to find a defect configuration
        infinite : bool (optional)
            Whether to create periodic infinite chain
        chain : Chain class (optional)
            Information about the polymer chain
        packing : str (optional)
            Chain packing method. Current support methods include: Hexagonal, Rectangular P, Rectangular C
        density : float (optional)
            Target packing density
        a : float (optional)
            Unit cell length
        b : float (optional)
            Unit cell length
        alpha : float (optional)
            Unit cell angle
        beta : float (optional)
            Unit cell angle
        gamma : float (optional)
            Unit cell angle
        """
        polymer_type = "PE"
        helice = Helice()
        num_monomers = 30
        tacticity = ""
        chiriality = ""
        head_tail_defect_ratio = 0
        configs = 30
        infinite = False

        for key in kwargs:
            if key == "polymer_type":
                polymer_type = kwargs["polymer_type"]
            elif key == "helice":
                helice = kwargs["helice"]
            elif key == "num_monomers":
                num_monomers = kwargs["num_monomers"]
            elif key == "tacticity":
                tacticity = kwargs["tacticity"]
            elif key == "chiriality":
                chiriality = kwargs["chiriality"]
            elif key == "head_tail_defect_ratio":
                head_tail_defect_ratio = kwargs["head_tail_defect_ratio"]
            elif key == "configs":
                configs = kwargs["configs"]
            elif key == "infinite":
                infinite = kwargs["infinite"]
            else:
                raise KeyError(
                    "Unknown input %s for Chain class\n Please see help for more information"
                    % key
                )

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
