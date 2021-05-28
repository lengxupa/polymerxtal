# ============================================================================
# scan.py -- Scanner functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import re

# from .stdio import FILE
from .utils import openFile

TOK_EOF = 0
TOK_INT = 1
TOK_REAL = 2
TOK_STRING = 3
TOK_DATA_DIR = 4
TOK_ELEMENT_DATA = 5
TOK_BOND_SCALE = 6
TOK_TEMPERATURE = 7
TOK_MONOMER = 8
TOK_TORSION = 9
TOK_FIXED = 10
TOK_FREE = 11
TOK_ENERGY = 12
TOK_STEREO = 13
TOK_PATTERN = 14
TOK_WEIGHT = 15
TOK_POLYMER = 16
TOK_BACKBONE_BOND_LENGTH = 17
TOK_DENSITY = 18
TOK_CHAINS = 19
TOK_MONOMERS = 20
TOK_MONOMERS_STDDEV = 21
TOK_SYSTEM_MAX = 22
TOK_SYSTEM_MIN = 23
TOK_EXCLUDE = 24
TOK_CYLINDER = 25
TOK_SLAB = 26
TOK_CHAIN_STEREO = 27
TOK_DOMAINS = 28
TOK_GRID_SIZE = 29
TOK_SAMPLE = 30
TOK_MONTE_CARLO = 31
TOK_NONE = 32
TOK_CONFIGS = 33
TOK_LJ = 34
TOK_ENERGY_CUTOFF = 35
TOK_BOND_CUTOFF = 36
TOK_RNG_SEED = 37
TOK_WRITE = 38
TOK_WRAPPED_PDB = 39
TOK_UNWRAPPED_PDB = 40
TOK_WRAPPED_XYZ = 41
TOK_UNWRAPPED_XYZ = 42
TOK_CHAIN_LENGTH_HISTO_BIN = 43
TOK_TORSION_ROTATION = 44
TOK_CHAIN_LENGTH = 45
TOK_CHAIN_LENGTH_HISTO = 46
TOK_TORSION_HISTO = 47
TOK_LOG_FILE = 48
TOK_STATUS_FILE = 49
TOK_EQUAL = 50
TOK_RECALCULATE_POSITIONS = 51
TOK_INTERMEDIATE = 52
TOK_RECALCULATE_NEIGHBORS = 53
TOK_INVERT = 54
TOK_TERM = 55
TOK_ISOLATE_CHAINS = 56
TOK_CALCULATE = 57
TOK_ALL = 58
TOK_TORSION_STEP = 59
TOK_SELF_AVOID = 60
TOK_SELF_AVOID_CUTOFF = 61
TOK_SPHERE = 62

switcher = {
    "=": TOK_EQUAL,
    "data_dir": TOK_DATA_DIR,
    "element_data ": TOK_ELEMENT_DATA,
    "bond_scale": TOK_BOND_SCALE,
    "temperature": TOK_TEMPERATURE,
    "monomer": TOK_MONOMER,
    "torsion": TOK_TORSION,
    "fixed": TOK_FIXED,
    "free": TOK_FREE,
    "energy": TOK_ENERGY,
    "stereo": TOK_STEREO,
    "pattern": TOK_PATTERN,
    "weight": TOK_WEIGHT,
    "polymer": TOK_POLYMER,
    "backbone_bond_length": TOK_BACKBONE_BOND_LENGTH,
    "density": TOK_DENSITY,
    "chains": TOK_CHAINS,
    "monomers": TOK_MONOMERS,
    "monomers_stddev": TOK_MONOMERS_STDDEV,
    "system_min": TOK_SYSTEM_MIN,
    "system_max": TOK_SYSTEM_MAX,
    "exclude": TOK_EXCLUDE,
    "cylinder": TOK_CYLINDER,
    "slab": TOK_SLAB,
    "sphere": TOK_SPHERE,
    "chain_stereo": TOK_CHAIN_STEREO,
    "domains": TOK_DOMAINS,
    "grid_size": TOK_GRID_SIZE,
    "sample": TOK_SAMPLE,
    "monte_carlo": TOK_MONTE_CARLO,
    "none": TOK_NONE,
    "configs": TOK_CONFIGS,
    "LJ": TOK_LJ,
    "energy_cutoff": TOK_ENERGY_CUTOFF,
    "bond_cutoff": TOK_BOND_CUTOFF,
    "rng_seed": TOK_RNG_SEED,
    "write": TOK_WRITE,
    "wrapped_pdb": TOK_WRAPPED_PDB,
    "unwrapped_pdb": TOK_UNWRAPPED_PDB,
    "wrapped_xyz": TOK_WRAPPED_XYZ,
    "unwrapped_xyz": TOK_UNWRAPPED_XYZ,
    "chain_length_histo_bin": TOK_CHAIN_LENGTH_HISTO_BIN,
    "torsion_rotation": TOK_TORSION_ROTATION,
    "chain_length": TOK_CHAIN_LENGTH,
    "chain_length_histo": TOK_CHAIN_LENGTH_HISTO,
    "torsion_histo": TOK_TORSION_HISTO,
    "intermediate": TOK_INTERMEDIATE,
    "log_file": TOK_LOG_FILE,
    "status_file": TOK_STATUS_FILE,
    "recalculate_positions": TOK_RECALCULATE_POSITIONS,
    "recalculate_neighbors": TOK_RECALCULATE_NEIGHBORS,
    "invert": TOK_INVERT,
    "term": TOK_TERM,
    "isolate_chains": TOK_ISOLATE_CHAINS,
    "calculate": TOK_CALCULATE,
    "all": TOK_ALL,
    "torsion_step": TOK_TORSION_STEP,
    "self_avoid": TOK_SELF_AVOID,
    "self_avoid_cutoff": TOK_SELF_AVOID_CUTOFF,
    "[-+]?[0-9]+": TOK_INT,
    "[-+]?\d*\.\d+|\d+": TOK_REAL,  # noqa: W605
    "[a-zA-Z0-9_/\\.-]+": TOK_STRING,
}


def yylex(item):
    scan_flag = 1
    for key in switcher:
        if re.match("^" + key + "$", item):
            tokval = switcher[key]
            tokstr = item
            scan_flag = 0
            break
    if scan_flag:
        return -1, "Invalid argument"
    return tokval, tokstr


class Scanner:
    def __init__(self):
        self.path = ""
        self.tokstr = ""
        self.scanner = ""
        self.tokval = 0
        self.lineno = 0
        self.reuse = 0
        self.toklen = 0

    # ============================================================================
    # initScanner()
    # ----------------------------------------------------------------------------
    # Result: initialize a Scanner -- open FILE, create internal scanner; call
    # choke() on error
    # ============================================================================
    def initScanner(self):
        self.f = openFile(self.path, "r")
        self.lineno = 1
        self.toklen = 0
        self.reuse = 0
        self.colno = 0

    # ============================================================================
    # resetScanner()
    # ----------------------------------------------------------------------------
    # Result: reset s to read from the beginning of the input file
    # ============================================================================
    def resetScanner(self):
        if self:
            self.f.fclose()
            self.scanner = ""
            self.tokstr = ""
            self.initScanner()

    # ============================================================================
    # pushToken()
    # ----------------------------------------------------------------------------
    # Result: push back the most recent token so that the next call to get*Token()
    # returns the same token
    # ============================================================================
    def pushToken(self):
        if self:
            self.reuse += 1

    # ============================================================================
    # getToken()
    # ----------------------------------------------------------------------------
    # Result: return the TOK_ value of the next input token; call choke() on error
    # and updates s->tokstr, s->toklen, s->lineno
    # ============================================================================
    def getToken(self):
        if self.reuse:
            self.reuse = 0
        else:
            self.colno += 1
            if self.colno >= self.toklen:
                self.line = self.f.readline()
                if self.line == "":
                    return TOK_EOF
                self.lineno += 1
                self.tok_list = self.line.split("#")[0].split()
                self.toklen = len(self.tok_list)
                while not self.toklen:
                    self.line = self.f.readline()
                    if self.line == "":
                        return TOK_EOF
                    self.lineno += 1
                    self.tok_list = self.line.split("#")[0].split()
                    self.toklen = len(self.tok_list)
                self.colno = 0
            self.scanner = self.tok_list[self.colno]
        self.tokval, self.tokstr = yylex(self.scanner)
        return self.tokval

    # ============================================================================
    # getIntToken()
    # ----------------------------------------------------------------------------
    # Result: return the next input token as an integer; calls choke() if the next
    # token can't be converted to an int
    # ============================================================================
    def getIntToken(self):
        self.tokval = self.getToken()
        if self.tokval == TOK_INT or self.tokval == TOK_REAL:
            r = int(float(self.tokstr))
        else:
            raise TypeError(
                'Expected an integer on line %d of file %s, but found "%s"'
                % (self.lineno, self.path, self.tokstr)
            )
        return r

    # ============================================================================
    # getRealToken()
    # ----------------------------------------------------------------------------
    # Result: return the next input token as a Real; calls choke() if the next
    # token can't be converted to a Real
    # ============================================================================
    def getRealToken(self):
        self.tokval = self.getToken()
        if self.tokval == TOK_INT or self.tokval == TOK_REAL:
            r = float(self.tokstr)
        else:
            raise TypeError(
                'Expected a Real on line %d of file %s, but found "%s"'
                % (self.lineno, self.path, self.tokstr)
            )
        return r


# ============================================================================
# createScanner()
# ----------------------------------------------------------------------------
# Result: return a pointer to a newly allocated, initialized Scanner; call
# choke() if allocation fails
# ============================================================================
def createScanner(path):
    s = Scanner()

    s.path = path
    s.tokstr = ""
    s.scanner = ""
    s.initScanner()
    return s
