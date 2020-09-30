"""
This module is for functions which build chain helix.
"""

import numpy as np
import os

from .helice import Helice
from .io import input_polymod, run_polymod, write_pdb
from .particle import Sphere, Cluster, Translator, Rotator
from .unit import chain_periodicity
from .visualize import ovito_view

# Get the location of the current module
current_location = os.path.dirname(__file__)


def calculate_rotation(vk):
    z = np.array([0, 0, 1])
    c = cross(chvec, z)


def build_chain(helice, tacticity='isotactic', chiriality='right', num_monomers=40):

    # Build the path to the sample files.
    in_path = os.path.join(current_location, 'data', 'polymod_input', 'sample_chain.txt')
    monomer_path = os.path.join(current_location, 'data', 'pdb', 'Propionitrile.pdb')

    # Build polymod input file
    input_polymod(in_path, monomer_path, helice, tacticity, chiriality, num_monomers, 'run_polymod.txt')

    # Run polymod
    h = run_polymod()

    # Correct chain orientation
    unit, unit_distance, vk = chain_periodicity(h)
    chain_cluster = Cluster()
    for atom in h.pos:
        chain_cluster.add_particle(Sphere(h.pos[atom]))

    z = np.array([0, 0, 1])
    axis = np.cross(vk, z)
    theta = np.arccos(np.dot(z, vk))
    rotator = Rotator()
    chain_cluster.move(rotator, axis=axis, theta=theta)

    for i in range(h.num_atoms):
        h.pos[i] = chain_cluster.particles[i].center

    helix_name = "helix_%s_%s" % (helice, '+' if chiriality == 'right' else '-')
    write_pdb("%s.pdb" % helix_name, h.el_names, h.pos)

    # View chain structure
    ovito_view('%s.pdb' % helix_name, "%s_%s.png" % (helix_name, 'Front'), 'Front')
    ovito_view('%s.pdb' % helix_name, "%s_%s.png" % (helix_name, 'Top'), 'Top')
