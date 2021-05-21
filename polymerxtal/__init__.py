"""
polymerXtal
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal
structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer,
tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given
atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers
where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline
polymer simulations.
"""

import os

# Create tmp folder
if not os.path.exists(".tmp"):
    os.mkdir(".tmp")
if not os.path.exists(".tmp/bonds"):
    os.mkdir(".tmp/bonds")
if not os.path.exists(".tmp/types"):
    os.mkdir(".tmp/types")

# Add imports here
from .data import sample_chain  # noqa: F401
from .crystal import Helice, Chain, Cell  # noqa: F401

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
