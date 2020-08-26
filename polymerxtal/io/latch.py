"""
Functions for manipulating Polymer Modeler files.
"""


def generate_latch_input(polymer_custom_input):
    polymer_custom_output = polymer_custom_input + '_latch.in'
    polymer_type_custom = {}
    src = open(polymer_custom_input)
    des = open(polymer_custom_output, 'w')
    des.write('#\n')
    des.write('\n')
    des.write("# You shouldn't need this line if you build PolymerModeler\n")
    des.write('data_dir latch #/home/ben/work/PolymerModeler/src/latch/data\n')
    des.write('\n')
    des.write('# Use this scale to extend VdW radii when searching for bonds\n')
    des.write('bond_scale 1.2 \n')
    des.write('\n')
    des.write('# Temperature (K)\n')
    des.write('temperature 300\n')
    des.write('\n')
    des.write('# Make bonds of this length (in Angstroms) along the backbone\n')
    des.write('backbone_bond_length 1.53\n')
    des.write('\n')
    des.write('# Atoms connected by fewer than this many bonds do NOT feel pair interactions\n')
    des.write('bond_cutoff 4\n')
    des.write('\n')
    num_monomer = 0
    num_chain = 0
    polymer_type_custom['num_monomer'] = 0
    polymer_type_custom['num_chain'] = 0
    for line in src.readlines():
        valid_line = line.split('#')[0]
        ln = valid_line.split()
        if len(ln) < 1:
            continue
        if ln[0] == 'monomer_type':
            if len(ln) > 1:
                polymer_type_custom['num_monomer'] = int(ln[1])
                if polymer_type_custom['num_monomer'] < 1:
                    print('Please specify a positive number of custom monomer types you wish to build')
                    return False
            else:
                print('Please specify a positive number of custom monomer types you wish to build')
                return False
        if ln[0] == 'monomer':
            if 'monomer' not in polymer_type_custom:
                polymer_type_custom['monomer'] = {}
            num_monomer += 1
            if num_monomer <= polymer_type_custom['num_monomer']:
                polymer_type_custom['monomer'][num_monomer] = {}
                if len(ln[1:]) == 3:
                    polymer_type_custom['monomer'][num_monomer]['info'] = ln[1:]
                    if not os.path.exists('../' + ln[1]):
                        print('Could not find monomer PDB File %s' % ln[1])
                        return False
                    elif int(ln[2]) == int(ln[3]):
                        print('Monomer tail atom should not be the same with the head atom')
                        return False
                    else:
                        des.write('#monomer: name PDB head tail\n')
                        des.write('monomer m%d ../%s %d %d\n' % (num_monomer, ln[1], int(ln[2]), int(ln[3])))
                        des.write('\n')
                        des.write('# Let all backbone torsions rotate when monomer m1 is used in a chain\n')
                        des.write('torsion all free \n')
                        des.write('\n')
                else:
                    print('Please specify monomer %d information with monomer PDB file, head atom & tail atom' %
                          num_monomer)
                    return False
            else:
                print('More monomer input than specified, Please change accordingly')
                return False
        if ln[0] == 'torsion':
            if 'torsion' not in polymer_type_custom['monomer'][num_monomer]:
                polymer_type_custom['monomer'][num_monomer]['torsion'] = {}
            if len(ln) > 2:
                polymer_type_custom['monomer'][num_monomer]['torsion'][ln[1]] = ln[2:]
                if int(ln[2]) < 0 or int(ln[2]) > 3:
                    print('Please specify backbone torsion angle probabilities type')
                    print('Starting with backbone atom 3, specify whether the torsion should change, and, if so, how.')
                    print('Values for specification:')
                    print('  0: no change')
                    print('  1: uniform probability for all angles')
                    print('  2: energy levels associated with specific angless')
                    print('  3: probability associated with specific angles')
                    return False
                elif int(ln[2]) == 2 or int(ln[2]) == 3:
                    if len(ln) > 3:
                        if not os.path.exists('../' + ln[3]):
                            print('Could not find file %s' % ln[3])
                            return False
                    elif int(ln[2]) == 2:
                        print('Please specify file for energy levels associated with specific angles')
                        return False
                    elif int(ln[2]) == 3:
                        print('Please specify file for probability associated with specific angles')
                        return False
            else:
                print('Please specify backbone torsion angle probabilities information')
                return False
        if ln[0] == 'chain_stereo':
            if len(ln) > 1:
                polymer_type_custom['num_chain'] = int(ln[1])
                if polymer_type_custom['num_chain'] < 1:
                    print('Please specify a positive number of custom polymer chain types you wish to build')
                    return False
                elif len(ln[2:]) == 2 * polymer_type_custom['num_chain']:
                    if 'chain' not in polymer_type_custom:
                        polymer_type_custom['chain'] = {}
                    for i in range(polymer_type_custom['num_chain']):
                        if int(ln[2 * i + 2]) not in polymer_type_custom['chain']:
                            polymer_type_custom['chain'][int(ln[2 * i + 2])] = {}
                        polymer_type_custom['chain'][int(ln[2 * i + 2])]['probability'] = eval(ln[2 * i + 3])
                else:
                    print('Please specify Distribution of %d chain types' % polymer_type_custom['num_chain'])
                    return False
            else:
                print('Please specify a positive number of custom polymer chain types you wish to build')
                return False
        if l[0] == 'chain_type':
            if 'chain' not in polymer_type_custom:
                polymer_type_custom['chain'] = {}
            num_chain += 1
            if len(l[1:]) > 1:
                if num_chain not in polymer_type_custom['chain']:
                    polymer_type_custom['chain'][num_chain] = {}
                if 'arrangement' not in polymer_type_custom['chain'][num_chain]:
                    polymer_type_custom['chain'][num_chain]['arrangement'] = {}
                if int(l[1]) == 0:
                    if int(l[2]) < 1:
                        print('Please specify a positive number of monomers in first pattern for chain type %d' %
                              num_chain)
                        return False
                    elif int(l[2]) != len(l[3:]):
                        print('Please specify a sequence of %d monomer(s) as repeated pattern for chain type %d' %
                              (int(l[2]), num_chain))
                        return False
                    else:
                        polymer_type_custom['chain'][num_chain]['arrangement'] = {
                            'type': 0,
                            'len': int(l[2]),
                            'sequence': l[3:]
                        }
                elif int(l[1]) == 1:
                    if polymer_type_custom['num_monomer'] != len(l[2:]):
                        print('Please specify %d probabilities of each monomer arrangement' % int(l[2]))
                        return False
                    else:
                        polymer_type_custom['chain'][num_chain]['arrangement'] = {'type': 1, 'sequence': l[2:]}
                else:
                    print('Please specify chain type %d information with arrangement: 0 = pattern, 1 = probability' %
                          num_chain)
                    return False
            else:
                print('Please specify chain type %d information with arrangement: 0 = pattern, 1 = probability' %
                      num_chain)
                return False
    src.close()
    des.write('# Pattern of monomers used to make chains\n')
    des.write('stereo s1 pattern 1 m1\n')
    des.write('\n')
    des.write('# System density in g/cm^3\n')
    des.write('density 0.5\n')
    des.write('\n')
    des.write('# Build this many chains\n')
    des.write('chains 0\n')
    des.write('\n')
    des.write('# Chains have this many monomers\n')
    des.write('monomers 0\n')
    des.write('\n')
    des.write('#exclude invert slab <real> <real> <real> <real> <real> <real>\n')
    des.write('\n')
    des.write('# Distribution of chain types: all chains of type s1\n')
    des.write('chain_stereo 1 s1 1.0\n')
    des.write('\n')
    des.write('# Size of system grid\n')
    des.write('grid_size 6.0\n')
    des.write('\n')
    des.write('# Chain packing algorithm\n')
    des.write('sample monte_carlo\n')
    des.write('\n')
    des.write('# Sample this many Monte Carlo configurations\n')
    des.write('configs 50\n')
    des.write('\n')
    des.write('# Energy expression: Lennard-Jones\n')
    des.write('energy LJ\n')
    des.write('\n')
    des.write('# Energy cutoff (Angstroms)\n')
    des.write('energy_cutoff 6.0\n')
    des.write('\n')
    des.write('# RNG seed; if commented out, or zero, use current time\n')
    des.write('#rng_seed 0\n')
    des.write('\n')
    des.write('# Output\n')
    des.write('# Write an output PDB file with NON-periodic atomic positions\n')
    des.write('#write unwrapped_pdb\n')
    des.write('#write wrapped_pdb | unwrapped_pdb | wrapped_xyz | unwrapped_xyz\n')
    des.write('\n')
    des.write("# Write intermediate files for input to Chunyu's code\n")
    des.write('write intermediate\n')
    des.write('\n')
    des.write('# Log output messages to this file; default is stdout\n')
    des.write('#log_file <path>\n')
    des.write('\n')
    des.write('# Log status messages to this file; default is stdout\n')
    des.write('#status_file <path>\n')
    des.write('\n')
    des.close()
    return polymer_custom_output, polymer_type_custom
