"""
Tests for io subpackage.
"""

import os, pytest
import numpy as np

import polymerxtal.io

# Get the location of the current module
current_location = os.path.dirname(__file__)


def test_read_lammps():
    """
    Test for reading a lammps input file.
    """

    # Build the path to the test file.
    test_file = os.path.join(current_location, '..', 'data', 'lmp_data', 'DIAMOND.data')

    # Use function to read lammps
    calculated_output = polymerxtal.io.read_data(test_file)

    # Put expected output based on input file.
    expected_output = {
        'atoms': 4.0,
        'atom types': 1.0,
        'xlo': 0.0,
        'xhi': 2.4534,
        'ylo': 0.0,
        'yhi': 4.2492,
        'zlo': 0.0,
        'zhi': 0.1,
        'Masses': {
            1: [12]
        },
        'Atoms': {
            1: [1, 0, 0.0, 0.0, 0.0],
            2: [1, 0, 2.4534, 1.4164, 0.0],
            3: [1, 0, 1.2267, 2.1246, 0.0],
            4: [1, 0, 1.2267, 3.541, 0.0]
        }
    }

    assert calculated_output == expected_output


def test_readbond_error():
    """
    Test for correct error return of readbond function
    """
    tmp_file = open('tmp_bonds.dat', 'w')
    tmp_file.write('1 1 1 2\n')
    tmp_file.close()

    with pytest.raises(Exception):
        polymerxtal.io.polymod.readbond('tmp_bonds.dat')

    os.remove('tmp_bonds.dat')


def test_validate_bonds_true():
    """
    Test for validate bonds are sucessfully connected
    """
    coords = [np.zeros(3), np.array([1, 0, 0])]

    tmp_file = open('tmp_bonds.dat', 'w')
    tmp_file.write('BONDS\n')
    tmp_file.write('\n')
    tmp_file.write('1 1 1 2\n')
    tmp_file.close()

    calculated_output = polymerxtal.io.polymod.validate_bonds(coords, 'tmp_bonds.dat')

    os.remove('tmp_bonds.dat')

    expected_output = True

    assert calculated_output == expected_output


def test_validate_bonds_false():
    """
    Test for validate bonds are not sucessfully connected
    """
    coords = [np.zeros(3), np.array([100, 0, 0])]

    tmp_file = open('tmp_bonds.dat', 'w')
    tmp_file.write('BONDS\n')
    tmp_file.write('\n')
    tmp_file.write('1 1 1 2\n')
    tmp_file.close()

    calculated_output = polymerxtal.io.polymod.validate_bonds(coords, 'tmp_bonds.dat')

    os.remove('tmp_bonds.dat')

    expected_output = False

    assert calculated_output == expected_output
