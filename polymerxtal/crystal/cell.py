"""
Functions to generate polymer crystals.
"""
import copy
import numpy as np

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
                Helice(2,1,1) - Planar zig-zag
                Helice(2,2,1)
                Helice(2,3,1)
                Helice(4,1,1)
                Helice(4,1,2) - Planar zig-zag (Syndiotactic)
                Helice(4,2,1)
                Helice(4,2,2)
                Helice(4,3,1)
                Helice(4,3,2)
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
            Chain packing method. Current support methods include: Hexagonal, RectangularP, RectangularC, Oblique, Square
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
        polymer_type = "PE"
        helice = Helice()
        num_monomers = 30
        tacticity = ""
        chiriality = ""
        head_tail_defect_ratio = 0
        configs = 30
        infinite = True
        packing = "Hexagonal"
        density = 1
        r_ab = np.sqrt(3)

        monomers_flag = 0
        packing_flag = 0
        density_flag = 0
        a_flag = 0
        b_flag = 0
        r_ab_flag = 0
        alpha_flag = 0
        beta_flag = 0
        gamma_flag = 0
        for key in kwargs:
            if key == "polymer_type":
                polymer_type = kwargs["polymer_type"]
            elif key == "helice":
                helice = kwargs["helice"]
            elif key == "num_monomers":
                num_monomers = kwargs["num_monomers"]
                monomers_flag = 1
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
            elif key == "chain":
                chain = kwargs["chain"]
            elif key == "packing":
                packing = kwargs["packing"]
                if packing not in [
                    "Oblique",
                    "RectangularP",
                    "RectangularC",
                    "Square",
                    "Hexagonal",
                ]:
                    raise Exception(
                        "Unknown packing method %s\nCurrent support packing methods include: Oblique, RectangularP (Rectangular Prism), RectangularC (Centered Rectangular), Square, Hexagonal"
                        % packing
                    )
                packing_flag = 1
            elif key == "density":
                density = kwargs["density"]
                density_flag = 1
            elif key == "a":
                a = kwargs["a"]
                a_flag = 1
            elif key == "b":
                b = kwargs["b"]
                b_flag = 1
            elif key == "r_ab":
                r_ab = kwargs["r_ab"]
                r_ab_flag = 1
            elif key == "alpha":
                alpha = kwargs["alpha"]
                alpha_flag = 1
            elif key == "beta":
                beta = kwargs["beta"]
                beta_flag = 1
            elif key == "gamma":
                gamma = kwargs["gamma"]
                gamma_flag = 1
            else:
                raise KeyError(
                    "Unknown input %s for Chain class\n Please see help for more information"
                    % key
                )

        if "chain" in kwargs:
            self.chain = chain
        else:
            if monomers_flag:
                self.chain = Chain(
                    polymer_type=polymer_type,
                    helice=helice,
                    num_monomers=num_monomers,
                    tacticity=tacticity,
                    chiriality=chiriality,
                    head_tail_defect_ratio=head_tail_defect_ratio,
                    configs=configs,
                    infinite=infinite,
                )
            else:
                self.chain = Chain(
                    polymer_type=polymer_type,
                    helice=helice,
                    tacticity=tacticity,
                    chiriality=chiriality,
                    head_tail_defect_ratio=head_tail_defect_ratio,
                    configs=configs,
                    infinite=infinite,
                )

        self.packing = packing
        self.density = density
        self.r_ab = r_ab

        if not self.chain.built:
            print("The Helix Chain is not configured. Start building chain ...")
            self.chain.build_chain()
        self.holder = self.chain.holder
        self.a = self.chain.x
        self.b = self.chain.y
        self.c = self.chain.z
        print(
            "Building chain finishes. Periodic length along the chain is " + str(self.c)
        )
        self.alpha = self.chain.alpha
        self.beta = self.chain.beta
        self.gamma = self.chain.gamma

        if a_flag:
            self.a = a
        if b_flag:
            self.b = b
        if a_flag and b_flag and r_ab_flag and self._r_ab != r_ab:
            raise ValueError(
                "Conflicting input cell parameters a and b with a/b ratio r_ab"
            )

        if alpha_flag:
            self.alpha = alpha
        if beta_flag:
            self.beta = beta
        if gamma_flag:
            self.gamma = gamma

        if packing_flag:

            # Check conflicts between Packing method and input cell parameters alpha & beta
            if self.alpha != 90:
                raise ValueError(
                    "Unit cell parameter alpha should be 90 for %s Packing"
                    % self.packing
                )
            elif self.beta != 90:
                raise ValueError(
                    "Unit cell parameter beta should be 90 for %s Packing"
                    % self.packing
                )

            # Check conflicts between Packing method and input cell parameters a, b or a/b ratio
            if self.packing == "RectangularP":
                if self.gamma != 90:
                    raise ValueError(
                        "Unit cell parameter gamma should be 90 for %s Packing"
                        % self.packing
                    )
            elif self.packing == "RectangularC":
                if self.gamma != 90:
                    if (a_flag and b_flag and self.a != self.b) or (
                        r_ab_flag and self.r_ab != 1
                    ):
                        raise ValueError(
                            "Unit cell parameters a and b should equal for %s Packing"
                            % self.packing
                        )
                    if not r_ab_flag:
                        self.r_ab = 1
            elif self.packing == "Square":
                if self.gamma != 90:
                    raise ValueError(
                        "Unit cell parameter gamma should be 90 for %s Packing"
                        % self.packing
                    )
                elif (a_flag and b_flag and self.a != self.b) or (
                    r_ab_flag and self.r_ab != 1
                ):
                    raise ValueError(
                        "Unit cell parameters a and b should equal for %s Packing"
                        % self.packing
                    )
                if not r_ab_flag:
                    self.r_ab = 1
            elif self.packing == "Hexagonal":
                if self.gamma == 90:
                    if (
                        a_flag
                        and b_flag
                        and (self._r_ab != np.sqrt(3) or self._r_ab != np.sqrt(1 / 3))
                    ) or (
                        r_ab_flag
                        and (self.r_ab != np.sqrt(3) or self.r_ab != np.sqrt(1 / 3))
                    ):
                        raise ValueError(
                            "Unit cell parameters a and b should equal to square root of 3 for %s Packing"
                            % self.packing
                        )
                elif self.gamma == 60 or self.gamma == 120:
                    if (a_flag and b_flag and self.a != self.b) or (
                        r_ab_flag and self.r_ab != 1
                    ):
                        raise ValueError(
                            "Unit cell parameters a and b should equal for %s Packing"
                            % self.packing
                        )
                    if not r_ab_flag:
                        self.r_ab = 1
                else:
                    raise ValueError(
                        "Unit cell parameter gamma should be 60, 90 or 120 for %s Packing"
                        % self.packing
                    )

        volume = calculate_volume(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        current_density = calculate_density(self.chain.holder.el_names, volume)
        if (
            self.packing == "RectangularC" or self.packing == "Hexagonal"
        ) and self.gamma == 90:
            current_density *= 2

        if a_flag and b_flag:
            if density_flag and self.density != current_density:
                raise ValueError(
                    "Input cell parameters do not match input target packing density %f g/cc"
                    % self.density
                )
        elif a_flag:
            self.b *= current_density / self.density
            if density_flag:
                if r_ab_flag and self.r_ab != self._r_ab:
                    raise ValueError(
                        "Input cell parameter a and a/b ratio does not match input target packing density %f g/cc"
                        % self.density
                    )
            else:
                self.b = self.a / self.r_ab
        elif b_flag:
            self.a *= current_density / self.density
            if density_flag:
                if r_ab_flag and self.r_ab != self._r_ab:
                    raise ValueError(
                        "Input cell parameter b and a/b ratio does not match input target packing density %f g/cc"
                        % self.density
                    )
            else:
                self.b = self.a / self.r_ab
        else:
            self.b = np.sqrt(
                self.a * self.b * current_density / self.density / self.r_ab
            )
            self.a = self.r_ab * self.b

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

    def build_unit_cell(self, **kwargs):
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
            self.packing == "RectangularC" or self.packing == "Hexagonal"
        ) and self.gamma == 90:
            holder = self.create_centered_chain()
            self.holder = combine_holders(self.holder, holder)

        write_pdb(self.crystal_name + ".pdb", self.holder.el_names, self.holder.pos)

        use_visualize = False
        create_lmpdata_file = False
        bondscale = 1.1
        ffield = "Dreiding"
        charge = "Gasteiger"
        create_lmpinput_file = False

        for key in kwargs:
            if key == "use_visualize":
                use_visualize = kwargs["use_visualize"]
            elif key == "create_lmpdata_file":
                create_lmpdata_file = kwargs["create_lmpdata_file"]
            elif key == "bondscale":
                bondscale = kwargs["bondscale"]
            elif key == "ffield":
                ffield = kwargs["ffield"]
            elif key == "charge":
                charge = kwargs["charge"]
            elif key == "create_lmpinput_file":
                create_lmpinput_file = kwargs["create_lmpinput_file"]
            else:
                raise KeyError(
                    "Unknown input %s for build_chain function\n Please see help for more information"
                    % key
                )

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
