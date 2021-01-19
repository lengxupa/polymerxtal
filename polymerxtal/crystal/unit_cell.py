import numpy as np
from math import sqrt

from polymerxtal.io import write_pdb, write_lmp_ifile
from polymerxtal.polymod import createHolder
from polymerxtal.struct2lammps import Create_Data_File
from polymerxtal.visualize import ovito_view

from .chain import Chain
from .helice import Helice
from .holder import combine_holders
from .measure import calculate_volume
from .molecule import calculate_density
from .monomer import PolymerType, polymer_types
from .move import Sphere, Cluster, Translator

class ObliqueUnitCell:
    def __init__(self, a, b, gamma, chain=None, helice=None):
        self.a = a
        self.b = b
        self.alpha = 90
        self.beta = 90
        self.gamma = gamma

        if chain is None:
            # Create default chain object
            self.chain = Chain()
        if not self.chain.built:
            print("The Helix Chain is not configured. Start building chain ...")
            self.chain.build_chain()
        self.holder = self.chain.holder
        self.c = self.chain.z
        
        if helice is None:
            # Create default helice object
            self.helice = Helice()

        self.packing = "Oblique"

    def __str__(self):
        return f"{self.packing} Unit Cell\n a:\t {self.a}\nb:\t {self.b}\n alpha:\t {self.alpha}\nbeta:\t {self.beta}\ngamma:\t {self.gamma}\n"


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
            chain_cluster.add_particle(Sphere(self.holder.pos[atom]))
            holder.el_names[atom] = self.holder.el_names[atom]

        translator = Translator()
        chain_cluster.move(translator, array=array)

        for i in range(holder.num_atoms):
            holder.pos[i] = chain_cluster.particles[i].center

        return holder

    def build_unit_cell(self, use_visualize=False, create_lmpdata_file=False, create_lmpinput_file=False, bondscale=1.1, ffield="Dreiding", charge="Gasteiger"):
        """polymerxtal.crystal.cell.Cell.build_unit_cell function builds a unit cell of a polymer crystal with specified information

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
        self.crystal_name = self.chain.helix_name + "_" + self.packing

        # Build unit cell
        if (
            self.packing == "rectangularc" or self.packing == "hexagonal"
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
            )

        # View chain structure
        if use_visualize:
            write_pdb(
                self.crystal_name + "_view.pdb",
                self.holder.el_names,
                self.holder.pos,
                connect=False,
            )
            if create_lmpdata_file:
                ovito_view(
                    self.crystal_name + ".data",
                    self.crystal_name + "_Front.png",
                    view="Front",
                )
                ovito_view(
                    self.crystal_name + ".data",
                    self.crystal_name + "_Top.png",
                    view="Top",
                )
            else:
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


class RectangularPUnitCell(ObliqueUnitCell):
    def __init__(self, a, b):
        super().__init__(a, b, 90)
        self.packing = "Rectangular Prism"

class RectangularCUnitCell(ObliqueUnitCell):
    def __init__(self, a):
        super().__init__(a, a, 90)
        self.packing = "Centered Rectangle"

class SquareUnitCell(ObliqueUnitCell):
    def __init__(self, a):
        super().__init__(a, a, 90)
        self.packing = "Square"

class HexagonalUnitCell(ObliqueUnitCell):
    def __init__(self, a, gamma):

        possible_angles = {
                            90: sqrt(1/3),
                            60: 1, 
                            120: 1,
                         }   


        if gamma not in possible_angles:
            raise ValueError(f"Invalid gamma angle for hexagonal unit cell. Input value {gamma}. Value should be 60, 90, or 120.")

        # Calculate the value of b based on input of a and gamma.
        b = a / possible_angles[gamma]

        super().__init__(a, b, gamma)
        self.packing = "Hexagonal"




def Cell(cell_type, a, b=None, gamma=None):
    """
    Factory method for returning a cell object.
    """

    # Check cell type and parameters
    cell_type = cell_type.lower()

    cells = {
        "oblique": [ObliqueUnitCell, [a, b, gamma]],
        "rectangularp": [RectangularPUnitCell, [a, b]],
        "rectangularc": [RectangularCUnitCell, [a]],
        "square": [SquareUnitCell, [a]],
        "hexagonal": [HexagonalUnitCell, [a, gamma]],
    }

    if cell_type not in cells:
        raise ValueError("Packing type {cell_type} not supported.")
    
    # Check angle restriction
    if cell_type  in ["rectangularp", "rectangularc", "square"] and gamma != 90:
        if gamma and gamma !=90:
            raise ValueError(f"Invalid gamma value given for cell type {cell_type}. Gamma must be 90 degrees")
    
    if cell_type == "hexagonal":
        possible_angles = {
                        90: sqrt(1/3),
                        60: 1, 
                        120: 1,
                        }   


        if gamma not in possible_angles:
            raise ValueError(f"Invalid gamma angle for hexagonal unit cell. Input value {gamma}. Value should be 60, 90, or 120.")
            
        # Calculate the value of b based on input of a and gamma.
        if b and b != a / possible_angles[gamma]:
            raise ValueError("Invalid value for b parameter in {cell_type} unit cell")  

    # Check cell size description for square unit cell
    if (cell_type == "square") and b:
        sides = [a, b]

        sides = set([x for x in sides if x])

        if len(sides) > 1:
            raise ValueError("Incorrect cell parameters for cell type {cell_type}")
    
    elif (cell_type == "rectangularc") and b:
        if a != b:
            raise ValueError(f"Invalid value for b parameter and cell type {cell_type}")

    
    cell_entry = cells[cell_type]

    return cell_entry[0](*cell_entry[1])


    



