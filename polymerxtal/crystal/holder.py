"""
Temporary structure used to hold atomic positions and atom types while creating and filling a crystal structure
"""

from polymerxtal.io import readbond
from polymerxtal.polymod import createHolder


def delete_atoms(holder, atom_list):
    for atom in atom_list:
        del holder.pos[atom - 1]
        del holder.el_names[atom - 1]
        del holder.atom_types[atom - 1]
        holder.num_atoms -= 1
    return holder


def reassign_atom_ids(holder, start_id=0):
    atom_ids = {}
    atom_id = 0
    holder_new = createHolder(holder.num_atoms)
    for i in sorted(holder.pos):
        holder_new.pos[atom_id] = holder.pos[i]
        holder_new.el_names[atom_id] = holder.el_names[i]
        holder_new.atom_types[atom_id] = holder.atom_types[i]
        atom_id += 1
        atom_ids[i + 1] = atom_id + start_id
    return holder_new, atom_ids


def build_bonds_data(
    holder,
    in_path=".tmp/bonds.dat",
    atom_remove_list=[],
    atom_link_list=[],
    out_path=".tmp/bonds.dat",
    sort_ids=False,
):
    holder = delete_atoms(holder, atom_remove_list)

    if sort_ids:
        holder, atom_ids = reassign_atom_ids(holder)

    bonds = readbond(in_path)

    bonds_file = open(out_path, "w")
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
        if sort_ids:
            bonds_file.write(
                "%d 1 %d %d\n" % (bond_id, atom_ids[atom1], atom_ids[atom2])
            )
        else:
            bonds_file.write("%d 1 %d %d\n" % (bond_id, atom1, atom2))
        bond_id += 1
    bonds_file.write(
        "%d 1 %d %d\n"
        % (bond_id, atom_ids[atom_link_list[0]], atom_ids[atom_link_list[1]])
    )
    bonds_file.close()
    return holder


# def delete_atoms(holder, atom_list):
#    for atom in atom_list:
#        del holder.pos[atom - 1]
#        del holder.el_names[atom - 1]
#        del holder.atom_types[atom - 1]
#        holder.num_atoms -= 1
#    return holder


def combine_holders(
    holder1,
    holder2,
    holder1_bond_path=".tmp/bonds.dat",
    holder2_bond_path=".tmp/bonds.dat",
    out_path=".tmp/bonds.dat",
):
    holder1, atom_ids1 = reassign_atom_ids(holder1)
    bonds1 = readbond(holder1_bond_path)
    holder2, atom_ids2 = reassign_atom_ids(holder2, start_id=holder1.num_atoms)
    bonds2 = readbond(holder2_bond_path)

    atom_id = 0
    holder_new = createHolder(holder1.num_atoms + holder2.num_atoms)
    for i in sorted(holder1.pos):
        holder_new.pos[atom_id] = holder1.pos[i]
        holder_new.el_names[atom_id] = holder1.el_names[i]
        holder_new.atom_types[atom_id] = holder1.atom_types[i]
        atom_id += 1
    for i in sorted(holder2.pos):
        holder_new.pos[atom_id] = holder2.pos[i]
        holder_new.el_names[atom_id] = holder2.el_names[i]
        holder_new.atom_types[atom_id] = holder2.atom_types[i]
        atom_id += 1

    bonds_file = open(out_path, "w")
    bonds_file.write("BONDS\n")
    bonds_file.write("\n")
    bond_id = 1
    for bond in bonds1:
        atom1 = bond[0]
        atom2 = bond[1]
        bonds_file.write("%d 1 %d %d\n" % (bond_id, atom_ids1[atom1], atom_ids1[atom2]))
        bond_id += 1
    for bond in bonds2:
        atom1 = bond[0]
        atom2 = bond[1]
        bonds_file.write("%d 1 %d %d\n" % (bond_id, atom_ids2[atom1], atom_ids2[atom2]))
        bond_id += 1
    bonds_file.close()

    return holder_new
