# ============================================================================
# params.py -- Params functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np
import os, random, time

from .config import MAX_BONDS, REAL_MAX
from .energy import setSelfAvoidCutoff, Energy
from .element import readElements
from .exclude import ExclCylinder, ExclSlab, ExclSphere
from .monomer import *
from .os import storeDir, changeDir, restoreDir
from .scan import *
from .stdio import FILE
from .stereo import Stereo, createStereo
from .types import writeAtomTypes, writeBondTypes
from .zmatrix import ZMatrix

# Get the location of the current module
current_location = os.path.dirname(__file__)

# Build the path to data_dir.
prog_datadir = os.path.join(current_location, "..", "data")

# File scope
read_elements = 0


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# Simulation parameters
class Params:
    # ============================================================================
    # initParams()
    # ----------------------------------------------------------------------------
    # Result: initilize p
    # ============================================================================
    def __init__(self):
        v0 = np.zeros(3)

        self.system_min = np.zeros(3)  # minimum point of simulation box
        self.system_max = np.zeros(3)  # maximum point of simulation box
        self.system_size = np.zeros(3)  # system_max - system_min
        self.domain_size = np.zeros(3)  # size of thread spatial domain
        self.known_monomers = Monomer()  # list
        self.known_stereo = Stereo()  # list
        self.chain_stereo = {}  # size num_stereo
        self.excluded_slabs = ExclSlab()  # list
        self.data_dir = prog_datadir  # directory name for data files
        self.element_data = "elements"  # path to element data file
        self.log_file = ""  # stdout by default
        self.status_file = ""  # stdout by default
        self.excluded_cylinders = ExclCylinder()  # list
        self.excluded_spheres = ExclSphere()  # list
        self.chain_stereo_weights = {}  # size num_stereo
        self.energy_func = Energy()  # energy expression
        self.bond_scale = 0.0  # scale equilibrium bond lengths to identify bonds
        self.temperature = 0.1  # K not zero
        self.density = 0.0  # g/cm^3
        self.grid_size = 0.0  # Angstroms
        self.energy_cutoff = 0.0  # Angstroms
        self.self_avoid_cutoff = 0.0  # Angstroms
        self.total_volume = 0.0  # available for packing chains; Angstroms^3
        self.excluded_volume = 0.0  # Angstroms^3
        self.torsion_step = 0.0  # degrees
        self.backbone_bond_length = 0.0  # Angstroms
        self.num_monomers_stddev = 0.0
        self.chain_length_histo_bin = 0.0  # Angstroms
        self.rng_seed = int(time.time())
        self.num_chains = 0
        self.num_stereo = 0  # size of chain_stereo[], chain_stereo_weights[]
        self.num_monomers = 0  # mean polymerization
        self.num_domains_x = 1  # spatial domains in X
        self.num_domains_y = 1  # spatial domains in Y
        self.num_domains_z = 1  # spatial domains in Z
        self.total_domains = 0
        self.num_configs = 30  # Sampling configurations
        self.bond_cutoff = 0  # max bonded interactions separation
        self.num_delta_steps = 0
        self.max_monomer_atoms = 0  # Most atoms in any Monomer
        self.recalculate_positions = 0  # flags
        self.recalculate_neighbors = 0
        self.sample_monte_carlo = 0
        self.explicit_volume = 0
        self.inverted_volume = 0
        self.isolate_chains = 0
        self.write_wrapped_pdb = 0
        self.write_unwrapped_pdb = 0
        self.write_wrapped_xyz = 0
        self.write_unwrapped_xyz = 0
        self.write_chain_length_histo = 0
        self.write_chain_length = 0
        self.write_torsion_histo = 0
        self.write_intermediate = 0
        self.half_system_size = np.zeros(3)

    # ============================================================================
    # getElements()
    # ----------------------------------------------------------------------------
    # Result: call readElements()
    # ============================================================================
    def getElements(self):
        global read_elements
        if not read_elements:
            leng = len(self.data_dir) + len(self.element_data) + 2

            el_path = "%s/%s" % (self.data_dir, self.element_data)
            readElements(el_path)
            el_path = ""
            read_elements += 1

    # ============================================================================
    # readParams()
    # ----------------------------------------------------------------------------
    # Result: fill p by reading indicated input file; calls choke() on input error
    # ============================================================================
    def readParams(self, path):
        s = createScanner(path)
        done = 0
        name = ""
        file_name = ""
        full_path = ""

        def CHOKE_PARSE(param):
            raise SyntaxError(
                "File %s, line %d: unknown %s, %s" % (s.path, s.lineno, param, s.tokstr)
            )

        while not done:
            tokval = s.getToken()
            if tokval == TOK_EOF:
                done += 1

            elif tokval == TOK_DATA_DIR:
                s.getToken()
                self.data_dir = ""
                self.data_dir = s.tokstr

            elif tokval == TOK_ELEMENT_DATA:
                s.getToken()
                self.element_data = ""
                self.element_data = s.tokstr

            elif tokval == TOK_BOND_SCALE:
                self.bond_scale = s.getRealToken()

            elif tokval == TOK_TEMPERATURE:
                self.temperature = s.getRealToken()

            elif tokval == TOK_BACKBONE_BOND_LENGTH:
                self.backbone_bond_length = s.getRealToken()

            elif tokval == TOK_DENSITY:
                self.density = s.getRealToken()

            elif tokval == TOK_CHAINS:
                self.num_chains = s.getIntToken()

            elif tokval == TOK_MONOMERS:
                self.num_monomers = s.getIntToken()

            elif tokval == TOK_MONOMERS_STDDEV:
                self.num_monomers_stddev = s.getRealToken()

            elif tokval == TOK_SYSTEM_MIN:
                self.system_min[0] = s.getRealToken()
                self.system_min[1] = s.getRealToken()
                self.system_min[2] = s.getRealToken()
                self.explicit_volume = 1

            elif tokval == TOK_SYSTEM_MAX:
                self.system_max[0] = s.getRealToken()
                self.system_max[1] = s.getRealToken()
                self.system_max[2] = s.getRealToken()
                self.explicit_volume = 1

            elif tokval == TOK_EXCLUDE:
                t = s.getToken()
                n = 0
                if TOK_INVERT == t:
                    n += 1
                    t = s.getToken()
                if t == TOK_CYLINDER:
                    ec = createExclCylinder()
                    ec.invert = n
                    if ec.invert:
                        self.inverted_volume = 1
                    ec.start[0] = s.getRealToken()
                    ec.start[1] = s.getRealToken()
                    ec.start[2] = s.getRealToken()
                    ec.axis[0] = s.getRealToken()
                    ec.axis[1] = s.getRealToken()
                    ec.axis[2] = s.getRealToken()
                    ec.radius = s.getRealToken()
                    ec.length = s.getRealToken()
                    ec.next = self.excluded_cylinders
                    self.excluded_cylinders = ec
                    self.excluded_volume += np.pi * ec.radius * ec.radius * ec.length
                elif t == TOK_SLAB:
                    es = createExclSlab()
                    es.invert = n
                    if es.invert:
                        self.inverted_volume = 1
                    es.min[0] = s.getRealToken()
                    es.min[1] = s.getRealToken()
                    es.min[2] = s.getRealToken()
                    es.max[0] = s.getRealToken()
                    es.max[1] = s.getRealToken()
                    es.max[2] = s.getRealToken()
                    if es.max[0] < es.min[0]:
                        raise ValueError("Excluded slab min x > max x")
                    if es.max[1] < es.min[1]:
                        raise ValueError("Excluded slab min y > max y")
                    if es.max[2] < es.min[2]:
                        raise ValueError("Excluded slab min z > max z")
                    es.next = self.excluded_slabs
                    self.excluded_slabs = es
                    self.excluded_volume += (
                        (es.max[0] - es.min[0])
                        * (es.max[1] - es.min[1])
                        * (es.max[2] - es.min[2])
                    )
                elif t == TOK_SPHERE:
                    esph = createExclSphere()
                    esph.invert = n
                    if esph.invert:
                        self.inverted_volume = 1
                    # TODO xcenter ycenter zcenter radius
                    esph.center[0] = s.getRealToken()
                    esph.center[1] = s.getRealToken()
                    esph.center[2] = s.getRealToken()
                    esph.radius = s.getRealToken()
                    esph.next = self.excluded_spheres
                    self.excluded_spheres = esph
                    self.excluded_volume += (
                        4.0 * np.pi * esph.radius * esph.radius * esph.radius / 3.0
                    )
                else:
                    CHOKE_PARSE("exclusion")

            elif tokval == TOK_GRID_SIZE:
                self.grid_size = s.getRealToken()

            elif tokval == TOK_DOMAINS:
                self.num_domains_x = s.getIntToken()
                self.num_domains_y = s.getIntToken()
                self.num_domains_z = s.getIntToken()

            elif tokval == TOK_LOG_FILE:
                s.getToken()
                self.log_file = s.tokstr

            elif tokval == TOK_STATUS_FILE:
                s.getToken()
                self.status_file = s.tokstr

            elif tokval == TOK_RNG_SEED:
                self.rng_seed = s.getIntToken()

            elif tokval == TOK_RECALCULATE_POSITIONS:
                self.recalculate_positions = 1

            elif tokval == TOK_RECALCULATE_NEIGHBORS:
                self.recalculate_neighbors = 1

            elif tokval == TOK_ENERGY_CUTOFF:
                self.energy_cutoff = s.getRealToken()

            elif tokval == TOK_SELF_AVOID_CUTOFF:
                self.self_avoid_cutoff = s.getRealToken()
                setSelfAvoidCutoff(self.self_avoid_cutoff)

            elif tokval == TOK_BOND_CUTOFF:
                self.bond_cutoff = s.getIntToken()
                if self.bond_cutoff > MAX_BONDS:
                    raise ValueError("bond_cutoff > MAX_BONDS")

            elif tokval == TOK_CONFIGS:
                self.num_configs = s.getIntToken()

            elif tokval == TOK_CHAIN_LENGTH_HISTO_BIN:
                self.chain_length_histo_bin = s.getRealToken()

            elif tokval == TOK_SAMPLE:
                t = s.getToken()
                if t == TOK_MONTE_CARLO:
                    self.sample_monte_carlo = 1
                elif t == TOK_NONE:
                    self.sample_monte_carlo = 0
                else:
                    CHOKE_PARSE("sample option")

            elif tokval == TOK_ENERGY:
                t = s.getToken()
                if t == TOK_LJ:
                    self.energy_func = Energy.energyLJ
                elif t == TOK_SELF_AVOID:
                    self.energy_func = Energy.energySelfAvoid
                elif t == TOK_NONE:
                    self.energy_func = Energy.energyNone
                else:
                    CHOKE_PARSE("energy option")

            elif tokval == TOK_WRITE:
                t = s.getToken()
                if t == TOK_WRAPPED_PDB:
                    self.write_wrapped_pdb = 1
                elif t == TOK_WRAPPED_XYZ:
                    self.write_wrapped_xyz = 1
                elif t == TOK_UNWRAPPED_PDB:
                    self.write_unwrapped_pdb = 1
                elif t == TOK_UNWRAPPED_XYZ:
                    self.write_unwrapped_xyz = 1
                elif t == TOK_TORSION_ROTATION:
                    s.getToken()
                    msearch = findMonomer(s.tokstr, self)
                    if (not msearch) or not hasattr(msearch, "next"):
                        CHOKE_PARSE("monomer")
                    a = s.getIntToken() - 1  # torsion index
                    b = s.getIntToken()  # angle step
                    leng = len(msearch.name) + 22
                    f = "%s_torsion_%2d_%3ddeg.pdb" % (msearch.name, a, b)
                    writeInternalRotationPDB(msearch, a, b, f)
                    f = ""
                elif t == TOK_CHAIN_LENGTH:
                    self.write_chain_length = 1
                elif t == TOK_CHAIN_LENGTH_HISTO:
                    self.write_chain_length_histo = 1
                elif t == TOK_TORSION_HISTO:
                    self.write_torsion_histo = 1
                elif t == TOK_INTERMEDIATE:
                    self.write_intermediate = 1
                else:
                    CHOKE_PARSE("write option")

            elif tokval == TOK_MONOMER:
                s.getToken()
                name = s.tokstr
                s.getToken()
                f = s.tokstr
                a = s.getIntToken() - 1  # head index
                b = s.getIntToken() - 1  # tail index
                self.getElements()
                m = readMonomer(name, f, a, b, self.bond_scale)
                m.next = self.known_monomers
                self.known_monomers = m
                name = ""
                f = ""
                if m.num_atoms > self.max_monomer_atoms:
                    self.max_monomer_atoms = m.num_atoms

            elif tokval == TOK_TORSION:
                if not (m and hasattr(m, "next")):
                    raise TypeError(
                        "File %s, line %d: no monomer specified for torsion"
                        % (s.path, s.lineno)
                    )
                if TOK_ALL == s.getToken():
                    a = 2
                    b = m.num_bb
                else:
                    s.pushToken()
                    b = s.getIntToken()
                    if b < 1 or b > m.num_bb:
                        raise ValueError(
                            "File %s, line %d: invalid torsion index %d"
                            % (s.path, s.lineno, b)
                        )
                    a = b - 1
                tmp = s.getToken()
                if tmp == TOK_FIXED:
                    for n in range(a, b):
                        m.torsions[n] = Torsion(0)
                        t = s.getToken()
                        if is_number(s.tokstr):
                            ang = float(s.tokstr)
                        else:  # not a number
                            s.pushToken()
                            ang = -REAL_MAX
                        for n in range(a, b):
                            m.setFixedTorsion(n, ang)
                elif tmp == TOK_FREE:
                    for n in range(a, b):
                        m.torsions[n] = Torsion(1)
                elif tmp == TOK_ENERGY:
                    for n in range(a, b):
                        m.torsions[n] = Torsion(2)
                    if TOK_CALCULATE == s.getToken():
                        s.getToken()
                        for n in range(a, b):
                            m.torsions[n] = TORSION_ENERGY_CALC
                            calculateTorsionEnergies(
                                m,
                                n,
                                self.bond_cutoff,
                                self.backbone_bond_length,
                                s.tokstr,
                            )
                            m.readTorsionEnergies(n, s.tokstr, self.temperature)
                    else:
                        for n in range(a, b):
                            m.readTorsionEnergies(n, s.tokstr, self.temperature)
                else:
                    CHOKE_PARSE("torsion option")

            elif tokval == TOK_STEREO:
                s.getToken()
                name = s.tokstr
                tmp = s.getToken()
                if tmp == TOK_PATTERN:
                    a = 1
                elif tmp == TOK_WEIGHT:
                    a = 0
                else:
                    CHOKE_PARSE("stereo option")
                b = s.getIntToken()  # number of monomers in the stereo
                st = createStereo(name, a, b)
                name = ""
                for n in range(b):
                    s.getToken()
                    msearch = findMonomer(s.tokstr, self)
                    if not (msearch and hasattr(msearch, "next")):
                        CHOKE_PARSE("monomer")
                    st.addStereoMonomer(msearch, 0.0 if a else s.getRealToken())
                if TOK_TERM == s.getToken():
                    s.getToken()
                    msearch = findMonomer(s.tokstr, self)
                    if not (msearch and hasattr(msearch, "next")):
                        CHOKE_PARSE("monomer")
                    st.term = msearch
                else:
                    s.pushToken()
                st.next = self.known_stereo
                self.known_stereo = st

            elif tokval == TOK_CHAIN_STEREO:
                self.num_stereo = s.getIntToken()
                for n in range(self.num_stereo):
                    s.getToken()
                    st = findStereo(s.tokstr, self)
                    if not (st and hasattr(st, "next")):
                        CHOKE_PARSE("stereo")
                    self.chain_stereo[n] = st
                    self.chain_stereo_weights[n] = s.getRealToken()

            elif tokval == TOK_POLYMER:
                s.getToken()
                self.getElements()
                q = s.tokstr.rfind("/")
                if q == -1:
                    raise ValueError(
                        "Expecting <polymer>/<torsion_option> polymer argument"
                    )
                leng = len(self.data_dir) + len(s.tokstr) + 11
                full_path = "%s/polymers/%s" % (self.data_dir, s.tokstr)
                storeDir()
                changeDir(full_path)
                self.readParams(s.tokstr[q + 1 :])
                restoreDir()
                full_path = ""

            elif tokval == TOK_ISOLATE_CHAINS:
                self.isolate_chains += 1

            elif tokval == TOK_TORSION_STEP:
                self.num_delta_steps = s.getIntToken()
                self.torsion_step = s.getRealToken()

            else:
                CHOKE_PARSE("keyword")
        del s

        # Sanity checks
        if self.temperature <= 0.0:
            raise ValueError("Invalid temperature: %f", self.temperature)
        if 0 == self.num_stereo:
            st = findStereo("default", self)
            if not (st and hasattr(st, "next")):
                raise TypeError("No chain stereo chemistry options specified")
            self.num_stereo = 1
            self.chain_stereo = {}
            self.chain_stereo[0] = st
            self.chain_stereo_weights = {}
            self.chain_stereo_weights[0] = 1.0

    # ============================================================================
    # reportParams()
    # ----------------------------------------------------------------------------
    # Result: write parameters to an output file
    # ============================================================================
    def reportParams(self, f):

        f.printf("Read data from %s\n" % self.data_dir)
        f.printf("Element info: %s/%s\n\n" % (self.data_dir, self.element_data))
        f.printf(
            "Scale equilibrium bond lengths by %f to identify bonds\n\n"
            % self.bond_scale
        )

        writeAtomTypes(f)
        writeBondTypes(f)
        f.printf("\n")

        m = self.known_monomers
        while m and hasattr(m, "next"):
            f.printf("Monomer %s:\n" % m.name)
            m.zm.writeZMatrix(f, 0)
            f.printf("Internal coordinates with Dreiding types:\n")
            m.zm.writeZMatrix(f, 1)
            f.printf("   %d backbone atoms\n" % m.num_bb)
            f.printf("   mass (without head and tail): %f amu\n" % m.central_mass)
            for i in range(2, m.num_bb):
                f.printf("   Backbone atom %d torsion angle: " % (i + 1))
                if m.torsions[i] == Torsion.TORSION_FIXED:
                    f.printf("fixed\n")
                elif m.torsions[i] == Torsion.TORSION_FREE:
                    f.printf("freely rotating\n")
                elif m.torsions[i] == Torsion.TORSION_ENERGY:
                    f.printf("bonded interactions (E(phi)) specified\n")
                elif m.torsions[i] == Torsion.TORSION_ENERGY_CALC:
                    f.printf("bonded interactions (E(phi)) calculated internally\n")
            f.printf(
                "   %d extra bonds not represented in z-matrix\n" % m.num_extra_bonds
            )
            b = m.extra_bonds
            while b and hasattr(b, "next"):
                f.printf(
                    "      Between atoms %d and %d\n" % (b.index1 + 1, b.index2 + 1)
                )
                b = b.next
            m = m.next

        f.printf("\nStereochemistry options:\n")
        s = self.known_stereo
        while s and hasattr(s, "next"):
            f.printf(
                "   %s (%s): "
                % (s.name, "pattern" if s.pattern else "weighted selection")
            )
            for i in range(s.num_monomers):
                f.printf("%s " % s.monomers[i].name)
                if not s.pattern:
                    f.printf("(%f) " % s.weights[i])
            if s.pattern:
                f.printf("(repeat)")
            f.printf("\n")
            s = s.next

        f.printf("\nSystem dimensions:\n")
        f.printf(
            "   Minimum (%f, %f, %f)\n"
            % (self.system_min[0], self.system_min[1], self.system_min[2])
        )
        f.printf(
            "   Maximum (%f, %f, %f)\n"
            % (self.system_max[0], self.system_max[1], self.system_max[2])
        )
        if (
            (self.excluded_cylinders and hasattr(self.excluded_cylinders, "next"))
            or (self.excluded_slabs and hasattr(self.excluded_slabs, "next"))
            or (self.excluded_spheres and hasattr(self.excluded_spheres, "next"))
        ):
            ec = self.excluded_cylinders
            es = self.excluded_slabs
            esph = self.excluded_spheres

            while ec and hasattr(ec, "next"):
                if ec.invert:
                    f.printf("   Excluded inverted")
                else:
                    f.printf("   Excluded")
                f.printf(
                    " cylinder: start at (%f, %f, %f), "
                    % (ec.start[0], ec.start[1], ec.start[2])
                )
                f.printf(
                    "axis (%f, %f, %f), radius %f, length %f\n"
                    % (ec.axis[0], ec.axis[1], ec.axis[2], ec.radius, ec.length)
                )
                ec = ec.next
            while es and hasattr(es, "next"):
                if es.invert:
                    f.printf("   Excluded inverted")
                else:
                    f.printf("   Excluded")
                    f.printf(
                        " slab: from (%f, %f, %f) to (%f, %f, %f)\n"
                        % (
                            es.min[0],
                            es.min[1],
                            es.min[2],
                            es.max[0],
                            es.max[1],
                            es.max[2],
                        )
                    )
                es = es.next
            while esph and hasattr(esph, "next"):
                if esph.invert:
                    f.printf("   Excluded inverted")
                else:
                    f.printf("   Excluded")
                f.printf(
                    " sphere: centered at (%f, %f, %f) with radius %g\n"
                    % (esph.center[0], esph.center[1], esph.center[2], esph.radius)
                )
                esph = esph.next
        f.printf("Total polymer volume: %g A^3\n\n" % self.total_volume)

        f.printf("Total domains: %d\n" % self.total_domains)
        f.printf(
            "   %d domain%s in X (size %f)\n"
            % (
                self.num_domains_x,
                "s" if self.num_domains_x > 1 else "",
                self.domain_size[0],
            )
        )
        f.printf(
            "   %d domain%s in Y (size %f)\n"
            % (
                self.num_domains_y,
                "s" if self.num_domains_y > 1 else "",
                self.domain_size[1],
            )
        )
        f.printf(
            "   %d domain%s in Z (size %f)\n"
            % (
                self.num_domains_z,
                "s" if self.num_domains_z > 1 else "",
                self.domain_size[2],
            )
        )

        f.printf("\nChain density: %f g/cm^3\n" % self.density)
        f.printf("%d chains\n" % self.num_chains)
        f.printf("%d monomers per chain\n" % self.num_monomers)
        f.printf("   %f deviation in monomers per chain\n\n" % self.num_monomers_stddev)
        for i in range(self.num_stereo):
            f.printf(
                "%f %% of chains will be %s\n"
                % (self.chain_stereo_weights[i] * 100.0, self.chain_stereo[i].name)
            )

        f.printf("\nChain growth:\n")
        f.printf("   consider %d configurations\n" % self.num_configs)
        if self.sample_monte_carlo:
            f.printf("   Monte Carlo sampling\n")
            f.printf("   Temperature: %f K\n" % self.temperature)
            f.printf("   Energy expression: ")
            if self.energy_func == Energy.energyLJ:
                f.printf("Lennard-Jones\n")
            else:
                f.printf("none\n")
            f.printf("   Bond cutoff for bonded interactions: %d\n" % self.bond_cutoff)
            f.printf(
                "   Interaction range for non-bonded interactions: %f A\n"
                % self.energy_cutoff
            )
            f.printf("   Grid size for neighbor bins: %f A\n" % self.grid_size)
        f.printf(
            "   Rotate varying torsions %f degrees when packing\n" % self.torsion_step
        )
        f.printf(
            "   Backbone bond length between monomers: %f A\n"
            % self.backbone_bond_length
        )

        f.printf("\nRNG seed: %d\n" % self.rng_seed)
        f.printf(
            "%s atomic positions\n"
            % ("Recalculate" if self.recalculate_positions else "Store")
        )
        if self.isolate_chains:
            f.printf("Grow isolated chains\n")
        if self.write_wrapped_pdb:
            f.printf("Write PDB file with folded positions\n")
        if self.write_unwrapped_pdb:
            f.printf("Write PDB file with unfolded positions\n")
        if self.write_wrapped_xyz:
            f.printf("Write XYZ file with folded positions\n")
        if self.write_unwrapped_xyz:
            f.printf("Write XYZ file with unfolded positions\n")
        if self.write_chain_length:
            f.printf("Write chain length\n")
        if self.write_chain_length_histo:
            f.printf(
                "Write chain length histogram (bin size %f A)\n"
                % self.chain_length_histo_bin
            )
        if self.write_torsion_histo:
            f.printf("Write torsion angle histogram\n")
        if self.write_intermediate:
            f.printf("Write intermediate files for LAMMPS preprocessing\n")


# ============================================================================
# findStereo()
# ----------------------------------------------------------------------------
# Result: return a pointer to the Stereo in the p->known_stereo list with
# matching name, else return NULL
# ============================================================================
def findStereo(name, p):
    st = p.known_stereo

    while st and hasattr(st, "next"):
        if name == st.name:
            break
        st = st.next
    return st


# ============================================================================
# findMonomer()
# ----------------------------------------------------------------------------
# Result: return a pointer to the Monomer in the p->known_monomers list with
# matching name, else return NULL
# ============================================================================
def findMonomer(name, p):
    m = p.known_monomers

    while m and hasattr(m, "next"):
        if name == m.name:
            break
        m = m.next
    return m
