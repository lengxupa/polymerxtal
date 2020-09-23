"""
This module is for functions which show sample helices.
"""

import os
from ovito.io import import_file
from ovito.vis import Viewport

from .helice import Helice

# Get the location of the current module
current_location = os.path.dirname(__file__)


def sample_chain(helice=Helice(), view='Perspective'):

    # Build the path to the sample file.
    sample_path = os.path.join(current_location, 'data', 'lmp_data',
                               'sample_helix_%d_%d_%d.data' % (helice.atoms, helice.motifs, helice.turns))

    # Import the sample file.
    pipeline = import_file(sample_path)
    pipeline.source.data.cell.vis.enabled = False
    pipeline.source.data.particles.vis.radius = .5
    pipeline.add_to_scene()

    vp = Viewport()

    if view == 'Perspective':
        vp.type = Viewport.Type.Perspective
        vp.camera_dir = (-1, -1, -1)
    elif view == 'Ortho':
        vp.type = Viewport.Type.Ortho
        vp.camera_dir = (-1, -1, -1)
    elif view == 'Top':
        vp.type = Viewport.Type.Top
    elif view == 'Bottom':
        vp.type = Viewport.Type.Bottom
    elif view == 'Front':
        vp.type = Viewport.Type.Front
    elif view == 'Back':
        vp.type = Viewport.Type.Back
    elif view == 'Left':
        vp.type = Viewport.Type.Left
    elif view == 'Right':
        vp.type = Viewport.Type.Right

    vp.zoom_all()
    vp.render_image(size=(800, 600),
                    filename="sample_helix_%d_%d_%d_%s.png" % (helice.atoms, helice.motifs, helice.turns, view),
                    background=(0, 0, 0),
                    frame=0)

    pipeline.remove_from_scene()
