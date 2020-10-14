"""
This module is for functions which build chain helix.
"""

import numpy as np
import os

from polymerxtal.data import polymer_types
from polymerxtal.io import run_polymod, write_pdb, readbond
from polymerxtal.visualize import ovito_view

from .helice import Helice
from .monomer import identify_backbone_path
from .move import Sphere, Cluster, Translator, Rotator
from .unit import chain_periodicity

# Get the location of the current module
current_location = os.path.dirname(__file__)


def calculate_rotation(vk):
    z = np.array([0, 0, 1])
    axis = np.cross(vk, z)
    theta = np.arccos(np.dot(z, vk))
    return axis, theta


def correct_chain_orientation(h, num_monomers):
    unit, unit_distance, vk = chain_periodicity(h, num_monomers)
    chain_cluster = Cluster()

    for atom in h.pos:
        chain_cluster.add_particle(Sphere(h.pos[atom]))

    axis, theta = calculate_rotation(vk)

    rotator = Rotator()
    chain_cluster.move(rotator, axis=axis, theta=theta)

    for i in range(h.num_atoms):
        h.pos[i] = chain_cluster.particles[i].center

    return h


def input_polymod(ofile,
                  monomers,
                  backbones,
                  tacticity='',
                  num_monomers=40):  #def input_polymod(nc, nm, cell, cpos, caxis, radius, lc, mctemp): origin

    des = open(ofile, 'w')
    des.write('\n')
    des.write('\n')
    des.write('bond_scale 1.2 \n')
    des.write('\n')
    des.write('temperature 300\n')
    des.write('\n')
    des.write('backbone_bond_length 1.5351\n')
    des.write('\n')
    des.write('bond_cutoff 4\n')
    des.write('\n')
    des.write('\n')

    for key in backbones:
        for T in ['60', '180', '300']:
            m = backbones[key].name + 'T' + T
            if m in monomers:
                des.write('monomer %s %s %d %d\n' % (backbones[key].name + 'T' + T, backbones[key].path,
                                                     backbones[key].head_index, backbones[key].tail_index))
                des.write('\n')
                des.write('torsion 3 fixed %s\n' % T)
                des.write('torsion 4 fixed %s\n' % backbones[key].name.split('T')[1])
                des.write('\n')
    des.write('\n')
    des.write('\n')
    #with open(in_path, 'r') as file:  #with open('file.txt', 'r') as file:
    #    filedata = file.read()

    # Replace the target string
    #filedata = filedata.replace('nc', str(nc))
    #filedata = filedata.replace('monomer_path', monomer_path)
    #filedata = filedata.replace('nm', str(num_monomers))
    #filedata = filedata.replace('lcx', str(cell[0]))
    #filedata = filedata.replace('lcy', str(cell[1]))
    #filedata = filedata.replace('lcz', str(cell[2]))
    #filedata = filedata.replace('ccx', str(cpos[0]))
    #filedata = filedata.replace('ccy', str(cpos[1]))
    #filedata = filedata.replace('ccz', str(cpos[2]))
    #filedata = filedata.replace('dx', str(caxis[0]))
    #filedata = filedata.replace('dy', str(caxis[1]))
    #filedata = filedata.replace('dz', str(caxis[2]))
    #filedata = filedata.replace('crad', str(radius))
    #filedata = filedata.replace('clc', str(lc + 200))
    #filedata = filedata.replace('irad', str(2 * radius - 5))
    #filedata = filedata.replace('mctemp', str(mctemp))

    stereo_type = ('%s 2 ' % ('weight' if tacticity == 'atactic' else 'pattern')) + monomers[0] + (
        ' 0.5 ' if tacticity == 'atactic' else ' ') + monomers[1] + (' 0.5' if tacticity == 'atactic' else '')
    #filedata = filedata.replace('stereo_type', stereo_type)
    des.write('stereo s1 %s\n' % stereo_type)
    des.write('\n')
    des.write('density 0.00001\n')
    des.write('\n')
    des.write('chains 1 #40\n')
    des.write('\n')
    des.write('\n')
    des.write('monomers %d\n' % num_monomers)
    des.write('\n')
    des.write('\n')
    des.write('\n')
    des.write('chain_stereo 1 s1 1.0\n')
    des.write('\n')
    des.write('grid_size 6.0\n')
    des.write('\n')
    des.write('sample monte_carlo\n')
    des.write('\n')
    des.write('configs 50\n')
    des.write('\n')
    des.write('energy LJ \n')
    des.write('\n')
    des.write('energy_cutoff 6.0\n')
    des.write('\n')
    des.write('\n')
    des.write('write unwrapped_pdb \n')
    des.write('\n')
    des.write('write intermediate\n')
    des.write('\n')
    des.write('\n')
    des.write('\n')
    des.write('\n')
    des.close()

    # Write the file out again
    #with open(ofile, 'w') as file:  #with open('run_file.txt', 'w') as file:
    #    file.write(filedata)


def find_configurations(monomer_path, helice, backbone_atoms=[], tacticity='', side_atom=0, chiriality=''):
    if len(backbone_atoms) != 2:
        raise ValueError('Number of backbone_atoms must be 2')
    if tacticity:
        if tacticity not in ['isotactic', 'atactic', 'syndiotactic']:
            raise TypeError(
                "Unknown tacticity, please specify one of the following: isotactic, atactic and syndiotactic")
        elif not side_atom:
            raise ValueError('Please specify side_atom')
    if side_atom and (not tacticity):
        raise ValueError('Please specify tacticity')

    if str(helice) == str(Helice(2, 3, 1)):
        if chiriality == 'left':
            torsions = ['T300', 'T180']
        elif chiriality == 'right':
            torsions = ['T60', 'T180']
        else:
            raise ValueError("Please specify chiriality: left or right")
    elif str(helice) == str(Helice(2, 2, 1)):
        if chiriality == 'left':
            torsions = ['T300']
        elif chiriality == 'right':
            torsions = ['T60']
        else:
            raise ValueError("Please specify chiriality: left or right")
    elif str(helice) == str(Helice(2, 1, 1)):
        torsions = ['T180']
        if chiriality:
            print('Zig-zag conformation does not have chiriality')

    backbones = identify_backbone_path(monomer_path, backbone_atoms, side_atom=side_atom)

    monomers = {}
    monomer_index = []

    for key in backbones:
        m = backbones[key].name.split('T')[0]
        if m not in monomers:
            monomers[m] = []
            monomer_index.append(m)
        for T in ['60', '180', '300']:
            m_name = backbones[key].name + 'T' + T
            check_flag = 1
            for tor in torsions:
                if tor not in m_name:
                    check_flag = 0
                    break
            if check_flag:
                monomers[m].append(m_name)

    configurations = {}
    config_index = 1
    if tacticity == 'isotactic' or (not tacticity):
        for m in monomer_index:
            for m_name in monomers[m]:
                configurations[config_index] = [m_name, m_name]
                config_index += 1
    elif tacticity == 'atactic' or tacticity == 'syndiotactic':
        for m1 in monomer_index:
            for m1_name in monomers[m1]:
                for m2 in monomer_index[monomer_index.index(m1) + 1:]:
                    for m2_name in monomers[m2]:
                        if m1_name[m1_name.index('T'):] == m2_name[m2_name.index('T'):]:
                            configurations[config_index] = [m1_name, m2_name]
                            config_index += 1

    return configurations, backbones


def validate_coords(coords, bond_path):
    bonds = readbond(bond_path)
    check_flag = 0
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            if ([i + 1, j + 1] not in bonds) and ([j + 1, i + 1] not in bonds):
                if np.linalg.norm(coords[i] - coords[j]) < 1.5:
                    check_flag = 1
                    return False
                    break
        if check_flag:
            break
    return True


def build_helix(monomer_path, helice, backbone_atoms=[], tacticity='', side_atom=0, chiriality='', num_monomers=40):

    configurations, backbones = find_configurations(monomer_path,
                                                    helice,
                                                    backbone_atoms=backbone_atoms,
                                                    tacticity=tacticity,
                                                    side_atom=side_atom,
                                                    chiriality=chiriality)

    print('Total %d possible configurations found for Helice %s' % (len(configurations), helice))

    for config_index in configurations:
        print('Try Configuration', config_index)

        # Build polymod input file
        input_polymod('run_polymod.txt',
                      configurations[config_index],
                      backbones,
                      tacticity=tacticity,
                      num_monomers=num_monomers)
        #input_polymod('run_polymod.txt',monomer_path,helice,backbone_atoms=backbone_atoms,tacticity=tacticity,side_atom=side_atom,chiriality=chiriality,num_monomers=num_monomers)
        #input_polymod('run_polymod.txt',monomer_path, helice, tacticity, chiriality, num_monomers=num_monomers)

        # Run polymod
        h = run_polymod('run_polymod.txt', validate_bond=True)

        if validate_coords(h.pos, 'bonds.dat'):
            print('Success for Configuration', config_index)
            break
        else:
            print('Fail for Configuration', config_index)
    return h


def build_chain(polymer_type, helice, num_monomers=40, tacticity='', chiriality='', use_visualize=False):

    # Build the path to the sample files.
    #in_path = os.path.join(current_location, '..', 'data', 'polymod_input', 'sample_chain.txt')
    #monomer_path = os.path.join(current_location, '..', 'data', 'pdb', 'monomer', 'PAN.pdb')

    if polymer_type not in polymer_types:
        raise ValueError(polymer_type + ' do not exist in our library, please consider using custom feature')

    # Build helix
    h = build_helix(polymer_types[polymer_type].path,
                    helice,
                    backbone_atoms=polymer_types[polymer_type].backbone_atoms,
                    tacticity=tacticity,
                    side_atom=polymer_types[polymer_type].side_atom,
                    chiriality=chiriality,
                    num_monomers=40)

    # Correct chain orientation
    h = correct_chain_orientation(h, num_monomers)

    helix_name = polymer_type + "_helix_%s_%s%s" % (helice, tacticity[:3], '+'
                                                    if chiriality == 'right' else '-' if chiriality == 'left' else '')
    write_pdb("%s.pdb" % helix_name, h.el_names, h.pos)

    # View chain structure
    if use_visualize:
        ovito_view('%s.pdb' % helix_name, "%s_%s.png" % (helix_name, 'Front'), 'Front')
        ovito_view('%s.pdb' % helix_name, "%s_%s.png" % (helix_name, 'Top'), 'Top')
