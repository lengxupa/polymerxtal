"""
unit.py
In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal
structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer,
tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given
atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers
where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline
polymer simulations.

Handles functions to calculate polymer crystals periodicity along the chain.
"""

import numpy as np


def normalized_dot_product(va, vb):
    return np.dot(va, vb) / np.linalg.norm(va) / np.linalg.norm(vb)


def chain_periodicity(h, num_monomers, unit_num_monomers):
    v0 = h.pos[3] - h.pos[1]

    num_atoms = int((h.num_atoms - 2) / num_monomers)

    indices = []
    dot_products = []
    distances = []
    for i in range(1, num_monomers):

        vi = h.pos[i * num_atoms + 3] - h.pos[(i - 1) * num_atoms + 3]

        indices.append((i - 1) * num_atoms + 3)
        dot_products.append(normalized_dot_product(v0, vi))
        distances.append(np.linalg.norm(h.pos[(i - 1) * num_atoms + 3] - h.pos[1]))

    distances_dir = {}

    for i, atom_id in enumerate(indices):
        if dot_products[i] > 0.99:
            distances_dir[i + 1] = distances[i]
        if distances[i] < 1:
            raise ValueError(
                "Distance of two monomers is less than 1 angstrom, indicating invalid chain structure"
            )

    k_flag = 0
    for k in sorted(distances_dir):
        if (
            (k * 2 in distances_dir)
            and (k * 3 in distances_dir)
            and (k * 5 in distances_dir)
        ):
            if -0.01 < (distances_dir[k * 2] - distances_dir[k] * 2) < 0.01:
                if -0.01 < (distances_dir[k * 3] - distances_dir[k] * 3) < 0.01:
                    if -0.01 < (distances_dir[k * 5] - distances_dir[k] * 5) < 0.01:
                        k_flag = 1
                        break

    if k_flag:
        unit = k
        distance = distances_dir[unit]
        max_unit = int(len(indices) / unit) * unit
        vk = h.pos[indices[max_unit - 1]] - h.pos[1]
    else:
        unit = unit_num_monomers
        # if distances_dir:
        #    unit = sorted(distances_dir)[0]
        # else:
        #    unit = num_monomers
        vk = h.pos[(unit - 1) * num_atoms + 3] - h.pos[1]
        distance = np.linalg.norm(vk)
    vk /= np.linalg.norm(vk)

    return unit, vk, distance
