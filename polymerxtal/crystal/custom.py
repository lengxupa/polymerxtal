import copy
import numpy as np

from polymerxtal.io import *  # noqa: F403

from polymerxtal.polymod import readPDB

def remove_atoms(holder,atoms):
	sorted_atoms=sorted(atoms)[::-1]
	num=1
	for atom_index in sorted_atoms:
		index=atom_index-num
		num+=1
		del holder.pos[index]
		del holder.el_names[index]
		del holder.atom_types[index]

		for atom in range(index,holder.num_atoms-1):
			holder.pos[atom]=copy.copy(holder.pos[atom+1])
			holder.el_names[atom]=copy.copy(holder.el_names[atom+1])
			holder.atom_types[atom]=copy.copy(holder.atom_types[atom+1])

		holder.num_atoms-=1
		del holder.pos[holder.num_atoms]
		del holder.el_names[holder.num_atoms]
		del holder.atom_types[holder.num_atoms]
		
	return holder

def add_holder(holder1,holder2,vector):
	holder=copy.copy(holder1)
	for atom in range(holder2.num_atoms):
		holder.pos[holder.num_atoms]=holder2.pos[atom]+vector
		holder.el_names[holder.num_atoms]=copy.copy(holder2.el_names[atom])
		holder.atom_types[holder.num_atoms]=copy.copy(holder2.atom_types[atom])
		holder.num_atoms+=1
	return holder


def custom_build_helix(polymer_type,num_monomers):
	unit_holder = readPDB(polymer_type.path)
	symbols, coords=open_pdb(polymer_type.path)

	key_atoms = polymer_type.backbone_atoms
	vector = unit_holder.pos[key_atoms[1] - 1] - unit_holder.pos[key_atoms[0] - 1]
	unit_distance = np.linalg.norm(vector)
	unit_vector = vector / unit_distance

	unit_head_atom = key_atoms[0]
	unit_head_append_atoms = polymer_type.append_atoms[unit_head_atom]
	unit_tail_atom = key_atoms[1]
	unit_tail_append_atoms = polymer_type.append_atoms[unit_tail_atom]

	for i in range(num_monomers):
		if i:			
			print('debug: custom.py: unit_holder', len(unit_holder.pos))
			print('debug: custom.py: holder', len(holder.pos))
			holder=remove_atoms(holder,tail_append_atoms)
			holder=add_holder(holder,unit_holder,i*vector)
			print('debug: custom.py: holder', len(holder.pos))
			tail_atom+=unit_holder.num_atoms
			for index in range(len(tail_append_atoms)):
				tail_append_atoms[index]+=1

		else:
			holder = copy.deepcopy(unit_holder)
			head_atoms=[unit_head_atom]+unit_head_append_atoms
			unit_holder=remove_atoms(unit_holder,head_atoms)
			print('debug: custom.py: unit_holder', unit_holder.num_atoms)
			print('debug: custom.py: holder', holder.num_atoms)
			tail_atom = unit_tail_atom
			tail_append_atoms = copy.copy(unit_tail_append_atoms)

	return holder, unit_vector, unit_distance

