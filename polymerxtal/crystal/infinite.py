"""
infinite.py
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal
structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer,
tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given
atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers
where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline
polymer simulations.

Handles functions to generate infinite periodic polymer crystal chain along chain direction
"""
import numpy as np

from .holder import build_bonds_data
from .move import Cluster, Sphere, Translator


def create_atom_list(holder, num_monomers):
    atom_remove_list = [1, 2, 3]

    num_atoms = int((holder.num_atoms - 2) / num_monomers)

    atom_remove_list += (
        [i for i in range(5, num_atoms + 3)]
        + [holder.num_atoms - 2 * num_atoms + 2]
        + [i for i in range(holder.num_atoms - num_atoms + 1, holder.num_atoms + 1)]
    )

    atom_link_list = [4, holder.num_atoms - 2 * num_atoms + 1]

    return atom_remove_list, atom_link_list


def create_infinite_chain(holder, num_monomers):
    atom_remove_list, atom_link_list = create_atom_list(holder, num_monomers)
    holder = build_bonds_data(
        holder,
        atom_remove_list=atom_remove_list,
        atom_link_list=atom_link_list,
        sort_ids=True,
    )
    return holder


def correct_infinite_chain_position(holder):
    trans_array = np.array([0, 0, -holder.pos[0][2]])
    chain_cluster = Cluster()

    for atom in holder.pos:
        chain_cluster.add_particle(Sphere(holder.pos[atom]))

    translator = Translator()
    chain_cluster.move(translator, array=trans_array)

    for i in range(holder.num_atoms):
        holder.pos[i] = chain_cluster.particles[i].center

    return holder
