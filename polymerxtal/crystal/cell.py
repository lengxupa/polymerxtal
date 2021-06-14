"""
Functions to generate polymer crystals.
"""
import copy
from math import sqrt
import numpy as np

from polymerxtal.io import write_pdb, write_lmp_ifile
from polymerxtal.polymod import createHolder
from polymerxtal.struct2lammps import Create_Data_File
from polymerxtal.visualize import ovito_view

from .chain import Chain

# from .helice import Helice
from .holder import combine_holders
from .measure import calculate_volume
from .molecule import calculate_density

# from .monomer import PolymerType, polymer_types
from .move import Sphere, Cluster, Translator


class UnitCell:
    def __init__(
        self,
        chain=None,
        packing=None,
        density=None,
        r_ab=None,
        a=None,
        b=None,
        alpha=None,
        beta=None,
        gamma=None,
    ):  # , helice=None):
        """polymerxtal.crystal.cell.UnitCell class contains all the information about the unit cell of polymer crystal

        Parameters
        ----------
        chain : Chain class (optional)
            Information about the polymer chain
        density : float (optional)
            Target packing density in grams per cubic centimeter
        r_ab : float (optional)
            a/b ratio
        a : float (optional)
            Unit cell length in Angstroms
        b : float (optional)
            Unit cell length in Angstroms
        alpha : float (optional)
            Unit cell angle in degrees
        beta : float (optional)
            Unit cell angle in degrees
        gamma : float (optional)
            Unit cell angle in degrees
        """
        # self.a = a
        # self.b = b
        # self.alpha = 90
        # self.beta = 90
        # self.gamma = gamma

        if chain is None:
            # Create default chain object
            self.chain = Chain(infinite=True)  # Chain()
        else:
            self.chain = chain

        if packing is None:
            self.packing = ""
        else:
            self.packing = packing
        if density is None:
            self.density = 1
        else:
            self.density = density
        if r_ab is None:
            self.r_ab = sqrt(3)
        else:
            self.r_ab = r_ab

        if not self.chain.built:
            print("The Helix Chain is not configured. Start building chain ...")
            self.chain.build_chain()
        self.holder = self.chain.holder

        if a is None:
            self.a = self.chain.x
        else:
            self.a = a
        if b is None:
            self.b = self.chain.y
        else:
            self.b = b
        if a and b and r_ab and self._r_ab != r_ab:
            raise ValueError(
                "Conflicting input cell parameters a and b with a/b ratio r_ab"
            )

        self.c = self.chain.z
        print(
            "Building chain finishes. Periodic length along the chain is " + str(self.c)
        )

        if alpha is None:
            self.alpha = self.chain.alpha
        else:
            self.alpha = alpha
        if beta is None:
            self.beta = self.chain.beta
        else:
            self.beta = beta
        if gamma is None:
            self.gamma = self.chain.gamma
        else:
            self.gamma = gamma

        # if helice is None:
        # Create default helice object
        #    self.helice = Helice()

        # self.packing = "Oblique"
        # self.packing = ""

        volume = calculate_volume(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        current_density = calculate_density(self.chain.holder.el_names, volume)
        if (
            self.packing == "Centered Rectangle" or self.packing == "Hexagonal"
        ) and self.gamma == 90:
            print("debug unitcell.py: UnitCell __init__: current_density")
            current_density *= 2

        if a and b:
            if density and self.density != current_density:
                raise ValueError(
                    "Input cell parameters do not match input target packing density %f g/cc"
                    % self.density
                )
        elif a:
            self.b *= current_density / self.density
            if density:
                if r_ab and self.r_ab != self._r_ab:
                    raise ValueError(
                        "Input cell parameter a and a/b ratio does not match input target packing density %f g/cc"
                        % self.density
                    )
            else:
                self.b = self.a / self.r_ab
        elif b:
            self.a *= current_density / self.density
            if density:
                if r_ab and self.r_ab != self._r_ab:
                    raise ValueError(
                        "Input cell parameter b and a/b ratio does not match input target packing density %f g/cc"
                        % self.density
                    )
            else:
                self.b = self.a / self.r_ab
        else:
            self.b = sqrt(self.a * self.b * current_density / self.density / self.r_ab)
            self.a = self.r_ab * self.b

    def __str__(self):
        return (
            "%s Unit Cell\n a:\t %f\nb:\t %f\nc:\t %f\n alpha:\t %f\nbeta:\t %f\ngamma:\t %f\n"
            % (self.packing, self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        )

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a):
        self._a = a
        if hasattr(self, "_b"):
            self._r_ab = self._a / self._b

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, b):
        self._b = b
        if hasattr(self, "_a"):
            self._r_ab = self._a / self._b

    def create_centered_chain(self):
        holder = createHolder(self.holder.num_atoms)

        array = np.array([self.a / 2, self.b / 2, 0])
        chain_cluster = Cluster()

        for atom in self.holder.pos:
            pos = copy.copy(self.holder.pos[atom])
            chain_cluster.add_particle(Sphere(pos))
            holder.el_names[atom] = self.holder.el_names[atom]

        translator = Translator()
        chain_cluster.move(translator, array=array)

        for i in range(holder.num_atoms):
            holder.pos[i] = chain_cluster.particles[i].center

        return holder

    def build_unit_cell(
        self,
        use_visualize=False,
        create_lmpdata_file=False,
        create_lmpinput_file=False,
        bondscale=1.1,
        ffield="Dreiding",
        charge="Gasteiger",
        lammps_min=False,
        lammps_min_levels=1,
    ):
        """polymerxtal.crystal.cell.UnitCell.build_unit_cell function builds a unit cell of a polymer crystal with
        specified information

        Parameters
        ----------
        use_visualize : bool (optional)
            Whether to view chain structure
        create_lmpdata_file : bool (optional)
            Whether to create corresponding LAMMPS data file or not
        bondscale : float (optional)
            Bond scale applied to equilibrium bond lengths
        ffield : str (optional)
            Force field. Current support force field include: Dreiding, PCFF
        charge : str (optional)
            Charge equilibration method. Current support charge method include: Gasteiger
        create_lmpinput_file : bool (optional)
            Whether to create corresponding LAMMPS input File or not

        """
        self.crystal_name = self.chain.helix_name + "_" + self.packing.replace(" ", "_")

        # Build unit cell
        if (
            self.packing == "Centered Rectangle" or self.packing == "Hexagonal"
        ) and self.gamma == 90:
            holder = self.create_centered_chain()
            self.holder = combine_holders(self.holder, holder)

        write_pdb(self.crystal_name + ".pdb", self.holder.el_names, self.holder.pos)

        # Create LAMMPS data file
        if create_lmpdata_file:
            Create_Data_File(
                self.crystal_name + ".pdb",
                bondscale=bondscale,
                ffield=ffield,
                charge=charge,
                xhi=self.a,
                yhi=self.b,
                zhi=self.c,
                alpha=self.alpha,
                beta=self.beta,
                gamma=self.gamma,
                outputName=self.crystal_name,
            )

            # Create LAMMPS input file
            if create_lmpinput_file:
                write_lmp_ifile(
                    ffield=ffield,
                    datafile=self.crystal_name + ".data",
                    potentialfile="X6paircoeffs.txt",
                    lammps_min=lammps_min,
                    lammps_min_levels=lammps_min_levels,
                )

        # View chain structure
        if use_visualize:
            write_pdb(
                self.crystal_name + "_view.pdb",
                self.holder.el_names,
                self.holder.pos,
                connect=False,
            )
            # if create_lmpdata_file:
            #    ovito_view(
            #        self.crystal_name + ".data",
            #        self.crystal_name + "_Front.png",
            #        view="Front",
            #    )
            #    ovito_view(
            #        self.crystal_name + ".data",
            #        self.crystal_name + "_Top.png",
            #        view="Top",
            #    )
            # else:
            # write_pdb(f"{helix_name}_ovito.pdb", h.el_names, h.pos, connect=False)
            ovito_view(
                self.crystal_name + "_view.pdb",
                self.crystal_name + "_Front.png",
                view="Front",
            )
            ovito_view(
                self.crystal_name + "_view.pdb",
                self.crystal_name + "_Top.png",
                view="Top",
            )

        return self.crystal_name


class ObliqueUnitCell(UnitCell):
    def __init__(self, chain=None, density=None, r_ab=None, a=None, b=None, gamma=None):
        super().__init__(
            chain=chain,
            packing="Oblique",
            density=density,
            r_ab=r_ab,
            a=a,
            b=b,
            alpha=90,
            beta=90,
            gamma=gamma,
        )
        # self.packing = "Oblique"


class RectangularPUnitCell(UnitCell):
    def __init__(self, chain=None, density=None, r_ab=None, a=None, b=None):
        super().__init__(
            chain=chain,
            packing="Rectangular Prism",
            density=density,
            r_ab=r_ab,
            a=a,
            b=b,
            alpha=90,
            beta=90,
            gamma=90,
        )
        # self.packing = "Rectangular Prism"


class RectangularCUnitCell(UnitCell):
    def __init__(self, chain=None, density=None, r_ab=None, a=None, b=None, gamma=None):
        if gamma and gamma != 90:
            # if a and b and a != b:
            # 	raise ValueError("Centered Rectangle Unit cell parameters a and b should equal")
            super().__init__(
                chain=chain,
                packing="Centered Rectangle",
                density=density,
                r_ab=1,
                a=a,
                b=b,
                alpha=90,
                beta=90,
                gamma=gamma,
            )
        else:
            super().__init__(
                chain=chain,
                packing="Centered Rectangle",
                density=density,
                r_ab=r_ab,
                a=a,
                b=b,
                alpha=90,
                beta=90,
                gamma=90,
            )
        # self.packing = "Centered Rectangle"


class SquareUnitCell(UnitCell):
    def __init__(self, chain=None, density=None, a=None, b=None):
        # if a and b and a != b:
        # 	raise ValueError("Square Unit cell parameters a and b should equal")
        super().__init__(
            chain=chain,
            packing="Square",
            density=density,
            r_ab=1,
            a=a,
            b=b,
            alpha=90,
            beta=90,
            gamma=90,
        )
        # self.packing = "Square"


class HexagonalUnitCell(UnitCell):
    def __init__(self, chain=None, density=None, r_ab=None, a=None, b=None, gamma=None):

        possible_angles = {90: [sqrt(3), sqrt(1 / 3)], 60: [1], 120: [1]}

        if gamma and (gamma not in possible_angles):
            raise ValueError(
                "Invalid gamma angle for hexagonal unit cell. Input value %f. Value should be 60, 90, or 120."
                % gamma
            )

        # elif gamma == 90:
        # 	if (a and b and (a/b != np.sqrt(3) or a/b != np.sqrt(1 / 3))) or (r_ab and (r_ab != np.sqrt(3) or r_ab !=
        #           np.sqrt(1 / 3))):
        # 		raise ValueError( "Hexagonal Unit cell parameters a/b ratio should equal to square root of 3 or 1/3")
        # 	super().__init__(chain = chain, density=density,r_ab=r_ab,a=a, b=b, alpha=90, beta=90, gamma=90)
        # else:
        # 	if a and b and a != b:
        # 		raise ValueError("Hexagonal Unit cell parameters a and b should equal")
        # 	super().__init__(chain = chain, density=density,r_ab=1,a=a, b=b, alpha=90, beta=90, gamma=gamma)

        # Calculate the value of b based on input of a and gamma.
        # b = a / possible_angles[gamma]

        # super().__init__(a, b, gamma)
        super().__init__(
            chain=chain,
            packing="Hexagonal",
            density=density,
            r_ab=r_ab,
            a=a,
            b=b,
            alpha=90,
            beta=90,
            gamma=gamma,
        )
        # self.packing = "Hexagonal"


def check_ratio(ratio, ratios):
    for r in ratios:
        if ratio == r:
            return True
    return False


def Cell(
    cell_type="hexagonal",
    chain=None,
    density=None,
    r_ab=None,
    a=None,
    b=None,
    alpha=None,
    beta=None,
    gamma=None,
):
    """
    Factory method for returning a cell object.
    """

    # Check cell type and parameters
    cell_type = cell_type.lower()

    cells = {
        "n/a": [UnitCell, [chain, density, r_ab, a, b, alpha, beta, gamma]],
        "oblique": [ObliqueUnitCell, [chain, density, r_ab, a, b, gamma]],
        "rectangularp": [RectangularPUnitCell, [chain, density, r_ab, a, b]],
        "rectangularc": [RectangularCUnitCell, [chain, density, r_ab, a, b, gamma]],
        "square": [SquareUnitCell, [chain, density, a, b]],
        "hexagonal": [HexagonalUnitCell, [chain, density, r_ab, a, b, gamma]],
    }

    if cell_type not in cells:
        raise ValueError("Packing type %s not supported." % cell_type)

    # Check angle restriction
    if cell_type in ["rectangularp", "square"] and gamma != 90:
        if gamma and gamma != 90:
            raise ValueError(
                "Invalid gamma value given for cell type %s. Gamma must be 90 degrees"
                % cell_type
            )

    if cell_type == "hexagonal":
        possible_angles = {
            90: [sqrt(3), sqrt(1 / 3)],
            60: [1],
            120: [1],
        }

        if gamma and (gamma not in possible_angles):
            raise ValueError(
                "Invalid gamma angle for hexagonal unit cell. Input value %f. Value should be 60, 90, or 120."
                % gamma
            )

        # Calculate the value of b based on input of a and gamma.
        if (a and b and (not check_ratio(a / b, possible_angles[gamma]))) or (
            r_ab and (not check_ratio(r_ab, possible_angles[gamma]))
        ):
            raise ValueError("Invalid value for a/b ratio in %s unit cell" % cell_type)

    # Check cell size description for square unit cell
    if (cell_type == "square") and a and b:
        sides = [a, b]

        sides = set([x for x in sides if x])

        if len(sides) > 1:
            raise ValueError(
                "Square Unit cell parameters a and b should equal"
            )  # "Incorrect cell parameters for cell type %s" %cell_type)

    elif (cell_type == "rectangularc") and a and b:
        if a != b and gamma and gamma != 90:
            raise ValueError(
                "Centered Rectangle Unit cell parameters a and b should equal"
            )  # "Incorrect cell parameters cell type %s" %cell_type)

    cell_entry = cells[cell_type]

    return cell_entry[0](*cell_entry[1])
