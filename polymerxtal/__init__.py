"""
polymerXtal
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer, tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline polymer simulations.
"""

# Add imports here
from .functions import *
from .io import open_pdb, open_xyz, write_xyz, write_lmp_ifile, run_lammps, read_cell_sizes, read_data, get_boundaries, read_input, check_args, complete_args
from .polymermodeler import main

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
