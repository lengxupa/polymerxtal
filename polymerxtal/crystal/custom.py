import copy
import numpy as np

from polymerxtal.io import *  # noqa: F403

from polymerxtal.polymod import readPDB

from .holder import build_bonds_data, combine_holders


def find_link_atoms(unit_holder, unit_head_atom, unit_head_append_atoms):
    bonds = readbond(".tmp/bonds.dat")  # noqa: F405
    link_atoms = []
    for bond in bonds:
        atom1 = bond[0]
        atom2 = bond[1]
        if atom1 == unit_head_atom and (atom2 not in unit_head_append_atoms):
            link_atoms.append(atom2)
        if atom2 == unit_head_atom and (atom1 not in unit_head_append_atoms):
            link_atoms.append(atom1)
    head_atoms = [unit_head_atom] + unit_head_append_atoms
    unit_holder = build_bonds_data(
        unit_holder,
        in_path=".tmp/bonds.dat",
        atom_remove_list=head_atoms,
        out_path=".tmp/unit_bonds.dat",
    )
    return unit_holder, link_atoms


def custom_build_helix(polymer_type, num_monomers, infinite):
    unit_holder = readPDB(polymer_type.path)
    symbols, coords = open_pdb(polymer_type.path)  # noqa: F405

    key_atoms = polymer_type.backbone_atoms
    vector = unit_holder.pos[key_atoms[1] - 1] - unit_holder.pos[key_atoms[0] - 1]
    unit_distance = np.linalg.norm(vector)
    unit_vector = vector / unit_distance

    unit_head_atom = key_atoms[0]
    unit_head_append_atoms = polymer_type.append_atoms[unit_head_atom]
    unit_tail_atom = key_atoms[1]
    unit_tail_append_atoms = polymer_type.append_atoms[unit_tail_atom]

    if infinite:
        for i in range(num_monomers):
            if i:
                holder = build_bonds_data(
                    holder,  # noqa: F405
                    in_path=".tmp/bonds.dat",
                    atom_remove_list=tail_append_atoms,  # noqa: F405
                    out_path=".tmp/bonds.dat",
                )
                atom_link_list = []
                for atom in link_atoms:  # noqa: F405
                    atom_link_list.append([tail_atom, atom])  # noqa: F405
                holder = combine_holders(
                    holder,
                    unit_holder,
                    holder1_bond_path=".tmp/bonds.dat",
                    holder2_bond_path=".tmp/unit_bonds.dat",
                    out_path=".tmp/bonds.dat",
                    vector=i * vector,
                    atom_link_list=atom_link_list,
                )
                tail_atom += unit_holder.num_atoms - num_tail_append_atoms  # noqa: F405
                for index in range(num_tail_append_atoms):  # noqa: F405
                    tail_append_atoms[index] += (  # noqa: F405
                        unit_holder.num_atoms - num_tail_append_atoms  # noqa: F405
                    )  # noqa: F405

            else:
                holder = copy.deepcopy(unit_holder)
                unit_holder, link_atoms = find_link_atoms(
                    unit_holder, unit_head_atom, unit_head_append_atoms
                )
                # head_atom = unit_head_atom
                head_append_atoms = copy.copy(unit_head_append_atoms)  # noqa: F841
                tail_atom = unit_tail_atom  # noqa: F841
                tail_append_atoms = copy.copy(unit_tail_append_atoms)  # noqa: F841
                num_tail_append_atoms = len(tail_append_atoms)

    else:
        for i in range(num_monomers):
            if i:
                holder = build_bonds_data(
                    holder,
                    in_path=".tmp/bonds.dat",
                    atom_remove_list=tail_append_atoms,
                    out_path=".tmp/bonds.dat",
                )
                atom_link_list = []
                for atom in link_atoms:
                    atom_link_list.append([tail_atom, atom])
                holder, atom_ids1, atom_ids2 = combine_holders(
                    holder,
                    unit_holder,
                    holder1_bond_path=".tmp/bonds.dat",
                    holder2_bond_path=".tmp/unit_bonds.dat",
                    out_path=".tmp/bonds.dat",
                    vector=i * vector,
                    atom_link_list=atom_link_list,
                    opt_ids=True,
                )
                tail_atom = atom_ids2[
                    unit_tail_atom
                ]  # -num_tail_append_atoms  # noqa: F405
                for index in range(num_tail_append_atoms):  # noqa: F405
                    tail_append_atoms[index] = atom_ids2[
                        unit_tail_append_atoms[index]
                    ]  # noqa: F405

            else:
                holder = copy.deepcopy(unit_holder)
                unit_holder, link_atoms = find_link_atoms(
                    unit_holder, unit_head_atom, unit_head_append_atoms
                )
                # head_atom = unit_head_atom
                head_append_atoms = copy.copy(unit_head_append_atoms)  # noqa: F841
                tail_atom = unit_tail_atom  # noqa: F841
                tail_append_atoms = copy.copy(unit_tail_append_atoms)  # noqa: F841
                num_tail_append_atoms = len(tail_append_atoms)

    return holder, unit_vector, unit_distance
