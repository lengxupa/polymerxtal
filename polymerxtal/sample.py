"""
This module is for functions which show sample helices.
"""

import os
from ovito.io import import_file
from ovito.vis import ParticlesVis, Viewport

from .helice import Helice

# Get the location of the current module
current_location = os.path.dirname(__file__)


def sample_chain(helice=Helice()):

    # Build the path to the sample file.
    sample_path = os.path.join(current_location, 'data', 'lmp_data',
                               'sample_helix_%d_%d_%d.data' % (helice.atoms, helice.motifs, helice.turns))

    # Import the sample file.
    pipeline = import_file(sample_path)
    pipeline.source.data.cell.vis.enabled = False
    pipeline.source.data.particles.vis.radius = .5
    pipeline.add_to_scene()

    vp = Viewport()
    vp.type = Viewport.Type.Front
    vp.zoom_all()
    vp.render_image(size=(800, 600),
                    filename="sample_helix_%d_%d_%d.png" % (helice.atoms, helice.motifs, helice.turns),
                    background=(0, 0, 0),
                    frame=0)

    pipeline.remove_from_scene()
