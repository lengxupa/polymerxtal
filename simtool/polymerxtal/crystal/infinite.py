"""
infinite.py
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer, tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline polymer simulations.

Handles functions to generate infinite periodic polymer crystal chain along chain direction
"""
import numpy as np
import os

from polymerxtal.io import readbond
from polymerxtal.polymod import createHolder

from .move import Cluster, Sphere, Translator


def create_atom_list(h, num_monomers):
    atom_remove_list = [1, 2, 3]

    num_atoms = int((h.num_atoms - 2) / num_monomers)

    atom_remove_list += (
        [i for i in range(5, num_atoms + 3)]
        + [h.num_atoms - 2 * num_atoms + 2]
        + [i for i in range(h.num_atoms - num_atoms + 1, h.num_atoms + 1)]
    )

    atom_link_list = [4, h.num_atoms - 2 * num_atoms + 1]

    return atom_remove_list, atom_link_list


def delete_atoms(h, atom_list):
    for atom in atom_list:
        del h.pos[atom - 1]
        del h.el_names[atom - 1]
        del h.atom_types[atom - 1]
        h.num_atoms -= 1
    return h


def reassign_atom_ids(h):
    atom_ids = {}
    atom_id = 0
    h_new = createHolder(h.num_atoms)
    for i in sorted(h.pos):
        h_new.pos[atom_id] = h.pos[i]
        h_new.el_names[atom_id] = h.el_names[i]
        h_new.atom_types[atom_id] = h.atom_types[i]
        atom_id += 1
        atom_ids[i + 1] = atom_id
    return h_new, atom_ids


def build_bonds_data(path, atom_remove_list, atom_link_list, atom_ids):
    bonds = readbond(path)

    if not os.path.exists(".tmp"):
        os.mkdir(".tmp")
    bonds_file = open(".tmp/bonds.dat", "w")
    bonds_file.write("BONDS\n")
    bonds_file.write("\n")
    bond_id = 1
    for bond in bonds:
        atom1 = bond[0]
        if atom1 in atom_remove_list:
            continue
        atom2 = bond[1]
        if atom2 in atom_remove_list:
            continue
        bonds_file.write("%d 1 %d %d\n" % (bond_id, atom_ids[atom1], atom_ids[atom2]))
        bond_id += 1
    bonds_file.write(
        "%d 1 %d %d\n"
        % (bond_id, atom_ids[atom_link_list[0]], atom_ids[atom_link_list[1]])
    )
    bonds_file.close()


def create_infinite_chain(h, num_monomers):
    atom_remove_list, atom_link_list = create_atom_list(h, num_monomers)
    h = delete_atoms(h, atom_remove_list)
    h, atom_ids = reassign_atom_ids(h)
    build_bonds_data(".tmp/bonds.dat", atom_remove_list, atom_link_list, atom_ids)
    return h


def correct_infinite_chain_position(h):
    trans_array = np.array([0, 0, -h.pos[0][2]])
    chain_cluster = Cluster()

    for atom in h.pos:
        chain_cluster.add_particle(Sphere(h.pos[atom]))

    translator = Translator()
    chain_cluster.move(translator, array=trans_array)

    for i in range(h.num_atoms):
        h.pos[i] = chain_cluster.particles[i].center

    return h
