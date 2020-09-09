# ============================================================================
# params.py -- Params functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import random
import time
import numpy as np

from .vector import Vector
from .monomer import Monomer
from .stereo import Stereo
from .exclude import ExclCylinder, ExclSlab, ExclSphere
from .energy import *
from .utils import FREE
from .scan import *
from .config import MAX_BONDS

# File scope
read_elements = 0


# ============================================================================
# findMonomer()
# ----------------------------------------------------------------------------
# Result: return a pointer to the Monomer in the p->known_monomers list with
# matching name, else return NULL
# ============================================================================
def findMonomer(name, p):
    m = p.known_monomers

    while m and hasattr(m, 'next'):
        if name == m.name:
            break
        m = m.next
    return m


# Simulation parameters
class Params:
    def __init__(self):
        self.system_min = Vector()  # minimum point of simulation box
        self.system_max = Vector()  # maximum point of simulation box
        self.system_size = Vector()  # system_max - system_min
        self.half_system_size = Vector()
        self.domain_size = Vector()  # size of thread spatial domain
        self.known_monomers = Monomer()  # list
        self.known_monomers.create()
        self.known_stereo = Stereo()  # list
        self.known_stereo.create()
        self.chain_stereo = {}  # size num_stereo
        self.excluded_cylinders = ExclCylinder()  # list
        self.excluded_cylinders.create()
        self.excluded_slabs = ExclSlab()  # list
        self.excluded_slabs.create()
        self.excluded_spheres = ExclSphere()  # list
        self.excluded_spheres.create()
        self.data_dir = ''  # directory name for data files
        self.element_data = "elements"  # path to element data file
        self.log_file = ''  # stdout by default
        self.status_file = ''  # stdout by default
        self.chain_stereo_weights = []  # size num_stereo
        self.bond_scale = 0.  # scale equilibrium bond lengths to identify bonds
        self.temperature = 0.1  # K not zero
        self.backbone_bond_length = 0.  # Angstroms
        self.density = 0.  # g/cm^3
        self.num_monomers_stddev = 0.
        self.grid_size = 0.  # Angstroms
        self.energy_cutoff = 0.  # Angstroms
        self.self_avoid_cutoff = 0.  # Angstroms
        self.chain_length_histo_bin = 0.  # Angstroms
        self.total_volume = 0.  # available for packing chains; Angstroms^3
        self.excluded_volume = 0.  # Angstroms^3
        self.torsion_step = 0.  # degrees
        self.energy_func = Energy()  # energy expression
        self.rng_seed = int(time.time())
        self.num_chains = 0
        self.num_stereo = 0  # size of chain_stereo[], chain_stereo_weights[]
        self.num_monomers = 0  # mean polymerization
        self.num_domains_x = 1  # spatial domains in X
        self.num_domains_y = 1  # spatial domains in Y
        self.num_domains_z = 1  # spatial domains in Z
        self.total_domains = 0
        self.max_monomer_atoms = 0  # Most atoms in any Monomer
        self.num_configs = 30  # Sampling configurations
        self.bond_cutoff = 0  # max bonded interactions separation
        self.num_delta_steps = 0
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

    # ============================================================================
    # getElements()
    # ----------------------------------------------------------------------------
    # Result: call readElements()
    # ============================================================================
    def getElements(self):
        if not read_elements:
            leng = len(self.data_dir) + len(self.element_data) + 2

            el_path = "%s/%s" % (self.data_dir, self.element_data)
            readElements(el_path)
            el_path = ''
            read_elements += 1

    # ============================================================================
    # freeParams()
    # ----------------------------------------------------------------------------
    # Result: free all allocated fields of p
    # ============================================================================
    def __del__(self):
        if self:
            es2 = ExclSlab()
            es2.create()

            FREE(self.chain_stereo_weights)
            self.data_dir = ''
            self.element_data = ''
            self.log_file = ''
            self.status_file = ''

            m1 = self.known_monomers
            while m1 and hasattr(m1, 'next'):
                m2 = m1.next
                del m1
                m1 = m2
            del self.known_monomers

            if self.chain_stereo:
                for i in range(self.num_stereo):
                    del self.chain_stereo[i]  # Free the Stereos below
                del self.chain_stereo
                self.chain_stereo = {}

            s1 = self.known_stereo
            while s1 and hasattr(s1, 'next'):
                s2 = s1.next
                del s1
                s1 = s2
            del self.known_stereo

            ec1 = self.excluded_cylinders
            while ec1 and hasattr(ec1, 'next'):
                ec2 = ec1.next
                del ec1
                ec1 = ec2
            del self.excluded_cylinders

            es1 = self.excluded_slabs
            while es1 and hasattr(es1, 'next'):
                es2 = es1.next
                del es1
                es1 = es2
            del self.excluded_slabs

    # ============================================================================
    # readParams()
    # ----------------------------------------------------------------------------
    # Result: fill p by reading indicated input file; calls choke() on input error
    # ============================================================================
    def readParams(self, path):
        s = Scanner().createScanner(path)
        done = 0
        ec = ExclCylinder()
        ec.create()
        es = ExclSlab()
        es.create()
        esph = ExclSphere()
        esph.create()
        msearch = Monomer()
        msearch.create()
        name = ''
        file_name = ''
        full_path = ''

        def CHOKE_PARSE(param):
            raise SyntaxError("File %s, line %d: unknown %s, %s" % (s.path, s.lineno, param, s.tokstr))

        while not done:
            tokval = s.getToken()
            if tokval == TOK_EOF:
                done += 1

            elif tokval == TOK_DATA_DIR:
                s.getToken()
                self.data_dir = ''
                self.data_dir = s.tokstr

            elif tokval == TOK_ELEMENT_DATA:
                s.getToken()
                self.element_data = ''
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
                self.system_min['x'] = s.getRealToken()
                self.system_min['y'] = s.getRealToken()
                self.system_min['z'] = s.getRealToken()
                self.explicit_volume = 1

            elif tokval == TOK_SYSTEM_MAX:
                self.system_max['x'] = s.getRealToken()
                self.system_max['y'] = s.getRealToken()
                self.system_max['z'] = s.getRealToken()
                self.explicit_volume = 1

            elif tokval == TOK_EXCLUDE:
                t = s.getToken()
                n = 0
                if TOK_INVERT == t:
                    n += 1
                    t = s.getToken()
                if t == TOK_CYLINDER:
                    ec = ExclCylinder()
                    ec.create()
                    ec.invert = n
                    if ec.invert:
                        self.inverted_volume = 1
                    ec.start.x = s.getRealToken()
                    ec.start.y = s.getRealToken()
                    ec.start.z = s.getRealToken()
                    ec.axis.x = s.getRealToken()
                    ec.axis.y = s.getRealToken()
                    ec.axis.z = s.getRealToken()
                    ec.radius = s.getRealToken()
                    ec.length = s.getRealToken()
                    ec.next = self.excluded_cylinders
                    self.excluded_cylinders = ec
                    self.excluded_volume += np.pi * ec.radius * ec.radius * ec.length
                elif t == TOK_SLAB:
                    es = ExclSlab()
                    es.create()
                    es.invert = n
                    if es.invert:
                        self.inverted_volume = 1
                    es.min.x = s.getRealToken()
                    es.min.y = s.getRealToken()
                    es.min.z = s.getRealToken()
                    es.max.x = s.getRealToken()
                    es.max.y = s.getRealToken()
                    es.max.z = s.getRealToken()
                    if es.max.x < es.min.x:
                        raise ValueError("Excluded slab min x > max x")
                    if es.max.y < es.min.y:
                        raise ValueError("Excluded slab min y > max y")
                    if es.max.z < es.min.z:
                        raise ValueError("Excluded slab min z > max z")
                    es.next = self.excluded_slabs
                    self.excluded_slabs = es
                    self.excluded_volume += (es.max.x - es.min.x) * (es.max.y - es.min.y) * (es.max.z - es.min.z)
                elif t == TOK_SPHERE:
                    esph = ExclSphere()
                    esph.create()
                    esph.invert = n
                    if esph.invert:
                        self.inverted_volume = 1
                    # TODO xcenter ycenter zcenter radius
                    esph.center.x = s.getRealToken()
                    esph.center.y = s.getRealToken()
                    esph.center.z = s.getRealToken()
                    esph.radius = s.getRealToken()
                    esph.next = self.excluded_spheres
                    self.excluded_spheres = esph
                    self.excluded_volume += 4.0 * np.pi * esph.radius * esph.radius * esph.radius / 3.0
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
                    self.energy_func = energyLJ
                elif t == TOK_SELF_AVOID:
                    self.energy_func = energySelfAvoid
                elif t == TOK_NONE:
                    self.energy_func = energyNone
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
                    if (not msearch) or not hasattr(msearch, 'next'):
                        CHOKE_PARSE("monomer")
                    a = s.getIntToken() - 1  # torsion index
                    b = s.getIntToken()  # angle step
                    leng = len(msearch.name) + 22
                    f = "%s_torsion_%2d_%3ddeg.pdb" % (msearch.name, a, b)
                    writeInternalRotationPDB(msearch, a, b, f)
                    f = ''
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
                m = readMonomer(name, f, a, b, self.bond_scale)  # Unfinished
                m.next = self.known_monomers
                self.known_monomers = m
                name = ''
                f = ''
                if m.num_atoms > self.max_monomer_atoms:
                    self.max_monomer_atoms = m.num_atoms

            #Unfinished

            elif tokval == TOK_TORSION:
                pass

            elif tokval == TOK_STEREO:
                pass

            elif tokval == TOK_CHAIN_STEREO:
                pass

            elif tokval == TOK_POLYMER:
                pass

            elif tokval == TOK_ISOLATE_CHAINS:
                pass

            elif tokval == TOK_TORSION_STEP:
                pass

            else:
                CHOKE_PARSE("keyword")


#void
#readParams(Params * restrict p, const char * restrict path)
#{
#   Scanner * restrict s = createScanner(path);
#   int done = 0;
#   int n, a, b, t;
#   ExclCylinder *ec;
#   ExclSlab *es;
#   ExclSphere *esph;
#   Stereo *st = NULL;
#   Monomer *m = NULL;
#   Monomer *msearch;
#   char *name = NULL;
#   char *file = NULL;
#   char *full_path = NULL;
#   char *q;
#   size_t len;
#   Real ang;
