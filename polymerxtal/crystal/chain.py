"""
This module is for functions Chain class.
"""

import numpy as np
import os, random
try:
    from ovito.io import import_file
    use_ovito = True
except ModuleNotFoundError:
    use_ovito = False

from polymerxtal.io import run_polymod, readbond, write_pdb
from polymerxtal.polymod import readPDB
from polymerxtal.visualize import ovito_view

from .helice import Helice
from .monomer import PolymerType, polymer_types
from .move import Sphere, Cluster, Translator, Rotator
from .unit import chain_periodicity


def validate_coords(coords_path, bond_path):
    h = readPDB(coords_path)
    bonds = readbond(bond_path)
    if use_ovito:
        pipeline = import_file(coords_path)
        types = pipeline.source.data.particles.particle_types
        for i in range(h.num_atoms):
            for j in range(i + 1, h.num_atoms):
                if ([i + 1, j + 1] not in bonds) and ([j + 1, i + 1] not in bonds):
                    if np.linalg.norm(h.pos[i] - h.pos[j]) <= (types.type_by_name(h.el_names[i]).radius +
                                                               types.type_by_name(h.el_names[j]).radius):
                        return False
    else:
        for i in range(h.num_atoms):
            for j in range(i + 1, h.num_atoms):
                if ([i + 1, j + 1] not in bonds) and ([j + 1, i + 1] not in bonds):
                    if np.linalg.norm(h.pos[i] - h.pos[j]) < 1.5:
                        return False
    return True


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


class Chain:
    def __init__(self,
                 polymer_type='PE',
                 helice=Helice(),
                 num_monomers=40,
                 tacticity='',
                 chiriality='',
                 head_tail_defect_ratio=0):

        if polymer_type not in polymer_types:
            raise ValueError(self.polymer_type + ' do not exist in our library, please consider using custom feature')

        if not 0 <= (head_tail_defect_ratio) <= 1:
            raise ValueError('Defect ratio of head to head and tail to tail connections is', head_tail_defect_ratio,
                             'and should be in the range of [0,1]')

        self.polymer_type = polymer_types[polymer_type]
        self.helice = helice
        if len(self.polymer_type.backbone_atoms) != self.helice.atoms:
            raise ValueError(f'Number of backbone_atoms must be {self.helice.atoms}')

        self.num_monomers = num_monomers
        self.tacticity = tacticity
        if self.tacticity:
            if self.tacticity not in ['isotactic', 'atactic', 'syndiotactic']:
                raise TypeError(
                    "Unknown tacticity, please specify one of the following: isotactic, atactic and syndiotactic")
            elif not self.polymer_type.backbone_atoms:
                raise ValueError('Please specify side_atom')

        self.chiriality = chiriality
        if str(self.helice) == str(Helice(2, 3, 1)):
            if self.chiriality == 'left':
                self.torsions = {180, 300}
            elif self.chiriality == 'right':
                self.torsions = {60, 180}
            else:
                raise ValueError("Please specify chiriality: left or right")
        elif str(self.helice) == str(Helice(2, 2, 1)):
            if self.chiriality == 'left':
                self.torsions = {300}
            elif self.chiriality == 'right':
                self.torsions = {60}
            else:
                raise ValueError("Please specify chiriality: left or right")
        elif str(self.helice) == str(Helice(2, 1, 1)):
            self.torsions = {180}
            if self.chiriality:
                print('Zig-zag conformation does not have chiriality')

        self.head_tail_defect_ratio = head_tail_defect_ratio
        self.pattern = 0
        self.monomers = []
        self.weights = {}

    def find_configurations(self, configs=50, find_inverse_path=False):

        if self.head_tail_defect_ratio:
            find_inverse_path = True
        self.backbones = self.polymer_type.identify_backbone_path(find_inverse_path=find_inverse_path)

        monomers = {}
        monomer_index = []

        for key in self.backbones:
            m = self.backbones[key].name.split('T')
            if m[0] not in monomers:
                monomers[m[0]] = []
                monomer_index.append(m[0])
            for T in [60, 180, 300]:
                m_name = self.backbones[key].name + 'T' + str(T)
                key_torsions = {eval(m[1]), T}
                if self.torsions == key_torsions:
                    monomers[m[0]].append(m_name)

        self.configurations = {}
        config_index = 1

        if find_inverse_path:
            self.pattern = 0
            while config_index <= configs:
                sequence = random.choices([0, 1], [1 - self.head_tail_defect_ratio, self.head_tail_defect_ratio],
                                          k=self.num_monomers - 1)
                if self.tacticity == 'isotactic' or (not self.tacticity):
                    m = random.choice(monomer_index)
                    m_name = random.choice(monomers[m])
                    h_t_label = m_name[0]
                    tor_label = m_name[m_name.index('T'):]
                    self.configurations[config_index] = [m_name]
                    for seq in sequence:
                        if seq:
                            m = random.choice(monomer_index)
                            m_name = random.choice(monomers[m])
                            while not (h_t_label != m_name[0] and tor_label == m_name[m_name.index('T'):]):
                                m = random.choice(monomer_index)
                                m_name = random.choice(monomers[m])
                            h_t_label = m_name[0]
                        self.configurations[config_index].append(m_name)
                elif self.tacticity == 'atactic' or self.tacticity == 'syndiotactic':
                    m1 = random.choice(monomer_index)
                    m1_name = random.choice(monomers[m1])
                    h_t_label = m1_name[0]
                    tor_label = m1_name[m1_name.index('T'):]
                    m2 = random.choice(monomer_index)
                    m2_name = random.choice(monomers[m2])
                    while not (h_t_label == m2_name[0] and m1 != m2 and tor_label == m2_name[m2_name.index('T'):]):
                        m2 = random.choice(monomer_index)
                        m2_name = random.choice(monomers[m2])
                    self.configurations[config_index] = [m1_name]
                    for i, seq in enumerate(sequence):
                        if seq:
                            m1 = random.choice(monomer_index)
                            m1_name = random.choice(monomers[m1])
                            while not (h_t_label != m1_name[0] and tor_label == m1_name[m1_name.index('T'):]):
                                m1 = random.choice(monomer_index)
                                m1_name = random.choice(monomers[m1])
                            h_t_label = m1_name[0]
                            m2 = random.choice(monomer_index)
                            m2_name = random.choice(monomers[m2])
                            while not (h_t_label == m2_name[0] and m1 != m2
                                       and tor_label == m2_name[m2_name.index('T'):]):
                                m2 = random.choice(monomer_index)
                                m2_name = random.choice(monomers[m2])
                        if i % 2:
                            self.configurations[config_index].append(m1_name)
                        else:
                            self.configurations[config_index].append(m2_name)
                config_index += 1
        else:
            if self.tacticity == 'isotactic' or (not self.tacticity):
                self.pattern = 0
                for m in monomer_index:
                    for m_name in monomers[m]:
                        self.configurations[config_index] = [m_name, m_name]
                        config_index += 1
            elif self.tacticity == 'atactic' or self.tacticity == 'syndiotactic':
                self.pattern = 1 if self.tacticity == 'atactic' else 0
                for m1 in monomer_index:
                    for m1_name in monomers[m1]:
                        for m2 in monomer_index[monomer_index.index(m1) + 1:]:
                            for m2_name in monomers[m2]:
                                if m1_name[m1_name.index('T'):] == m2_name[m2_name.index('T'):]:
                                    self.configurations[config_index] = [m1_name, m2_name]
                                    config_index += 1

    def input_polymod(self, ofile, config_index):

        self.monomers = set(self.configurations[config_index])

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

        for key in self.backbones:
            for T in [60, 180, 300]:
                m = self.backbones[key].name + 'T' + str(T)
                if m in self.monomers:
                    des.write(
                        f'monomer {m} {self.backbones[key].path} {self.backbones[key].head_index} {self.backbones[key].tail_index}\n'
                    )
                    des.write('\n')
                    des.write(f'torsion 3 fixed {T}\n')
                    des.write('torsion 4 fixed %s\n' % self.backbones[key].name.split('T')[1])
                    des.write('\n')
        des.write('\n')
        des.write('\n')
        #with open(in_path, 'r') as file:  #with open('file.txt', 'r') as file:
        #    filedata = file.read()

        # Replace the target string
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

        if self.pattern:
            weight = 1 / len(self.monomers)
            for m in self.monomers:
                self.weights[m] = weight
            stereo_type = 'weight ' + str(len(self.monomers))
        else:
            stereo_type = 'pattern ' + str(len(self.configurations[config_index]))

        for m in self.configurations[config_index]:
            stereo_type += ' ' + m
            if self.pattern:
                stereo_type += ' ' + str(self.weights[m])
        #filedata = filedata.replace('stereo_type', stereo_type)
        des.write('stereo s1 %s\n' % stereo_type)
        des.write('\n')
        des.write('density 0.00001\n')
        des.write('\n')
        des.write('chains 1 #40\n')
        des.write('\n')
        des.write('\n')
        des.write('monomers %d\n' % self.num_monomers)
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

    def build_helix(self):
        self.find_configurations()

        print('Total %d possible configurations found for Helice %s' % (len(self.configurations), self.helice))

        config_success = 0
        for config_index in self.configurations:
            print('Try Configuration', config_index)

            # Build polymod input file
            self.input_polymod('run_polymod.txt', config_index)
            #input_polymod('run_polymod.txt',monomer_path,helice,backbone_atoms=backbone_atoms,tacticity=tacticity,side_atom=side_atom,chiriality=chiriality,num_monomers=num_monomers)
            #input_polymod('run_polymod.txt',monomer_path, helice, tacticity, chiriality, num_monomers=num_monomers)

            # Run polymod
            run_polymod('run_polymod.txt', validate_bond=True)
            if validate_coords('chains_unwrapped.pdb', 'bonds.dat'):
                config_success = 1
                print('Success for Configuration', config_index)
                break
            else:
                print('Fail for Configuration', config_index)

        if config_success:
            return readPDB('chains_unwrapped.pdb')
        else:
            raise Exception(
                f'Unable to generate a successful Helice_{self.helice} configuration with {self.polymer_type.name} due to unavoided overlap of atoms.'
            )

    def build_chain(self, use_visualize=False):
        # Build the path to the sample files.
        #in_path = os.path.join(current_location, '..', 'data', 'polymod_input', 'sample_chain.txt')
        #monomer_path = os.path.join(current_location, '..', 'data', 'pdb', 'monomer', 'PAN.pdb')

        # Build helix
        h = self.build_helix()

        # Correct chain orientation
        h = correct_chain_orientation(h, self.num_monomers)

        helix_name = self.polymer_type.name + "_helix_%s%s%s%s" % (
            self.helice, '_' + self.tacticity[:3] if self.tacticity else '', '+' if self.chiriality == 'right' else
            '-' if self.chiriality == 'left' else '', '_custom' if self.head_tail_defect_ratio else '')
        write_pdb(f"{helix_name}.pdb", h.el_names, h.pos)

        # View chain structure
        if use_visualize:
            ovito_view(f'{helix_name}.pdb', f"{helix_name}_Front.png", 'Front')
            ovito_view(f'{helix_name}.pdb', f"{helix_name}_Top.png", 'Top')
