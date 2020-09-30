"""
This module is for functions which show sample helices.
"""

import os

from .helice import Helice
from .visualize import ovito_view

# Get the location of the current module
current_location = os.path.dirname(__file__)


def sample_chain(helice=Helice(), view='Perspective'):

    # Build the path to the sample file.
    sample_path = os.path.join(current_location, 'data', 'lmp_data',
                               'sample_helix_%d_%d_%d.data' % (helice.atoms, helice.motifs, helice.turns))

    ovito_view(sample_path, "sample_helix_%d_%d_%d_%s.png" % (helice.atoms, helice.motifs, helice.turns, view), view)

