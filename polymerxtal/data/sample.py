"""
This module is for functions which show sample helices.
"""

import os

from polymerxtal.crystal import Helice
from polymerxtal.visualize import ovito_view

# Get the location of the current module
current_location = os.path.dirname(__file__)


def sample_chain(helice=Helice(), view='Perspective'):

    # Build the path to the sample file.
    sample_path = os.path.join(current_location, 'lmp_data',
                               'sample_helix_%d_%d_%d.data' % (helice.atoms, helice.motifs, helice.turns))

    ovito_view(sample_path, "sample_helix_%d_%d_%d_%s.png" % (helice.atoms, helice.motifs, helice.turns, view), view)


class PolymerType:
    def __init__(self, name, path, backbone_atoms=[], side_atom=0):
        self.name = name
        self.path = path
        self.backbone_atoms = backbone_atoms
        self.side_atom = side_atom


PAN_path = os.path.join(current_location, 'pdb', 'monomer', 'PAN.pdb')
PE_path = os.path.join(current_location, 'pdb', 'monomer', 'PE.pdb')
PP_path = os.path.join(current_location, 'pdb', 'monomer', 'PP.pdb')

polymer_types = {
    'PAN': PolymerType('PAN', PAN_path, backbone_atoms=[3, 2], side_atom=5),
    'PE': PolymerType('PE', PE_path, backbone_atoms=[2, 3]),
    'PP': PolymerType('PP', PP_path, backbone_atoms=[2, 3], side_atom=8)
}
